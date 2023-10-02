import pybamm
import numpy as np
from math import sqrt
from matplotlib import pyplot as plt
import consts as c

class SingleParticle:

    def __init__(self, name: str, charge: int):
        self.name = name
        self.charge = charge
        self.domain = name + " dDomain"

        self.conc_0 = pybamm.Parameter(name + " pInitial Concentration")

        self.L = pybamm.Parameter(name + " pElectrode Thickness")
        self.eps_n = pybamm.Parameter(name + " pElectrode Porosity")

        self.conc_name = name + " vConcentration"
        self.flux_name = name + " vFlux"
        self.surf_conc_name = name + " vSurface Concentration"

        self.conc = pybamm.Variable(self.conc_name, domain=self.domain)
        self.surf_conc = pybamm.surf(self.conc)
        self.r = pybamm.SpatialVariable(name + " svRadius", domain=self.domain, coord_sys="spherical polar")

    def process_model(self, model: pybamm.BaseModel, current, clear=False):
        if clear:
            model = pybamm.BaseModel()

        flux = c.D * -pybamm.grad(self.conc)
        dcdt = -pybamm.div(flux)

        model.rhs.update({
            self.conc: dcdt,
        })

        model.initial_conditions.update({
            self.conc: self.conc_0,
        }) 
        
        a_term = (3 * (1 - self.eps_n)) / c.R

        self.j0 = pybamm.sqrt(self.surf_conc) * pybamm.sqrt(1 - self.surf_conc)
        self.j = current / (self.L * a_term)

        self.j_name = self.name + " vCurrent Density"

        model.boundary_conditions.update({
            self.conc: {
                "left":  (0, "Neumann"),
                "right": (self.j / (self.charge * c.F * c.D), "Neumann")
            },
        })
        model.variables.update({
            self.conc_name: self.conc, 
            self.surf_conc_name: self.surf_conc,
            self.j_name: self.j
        })
    
    def process_geometry(self, geo: dict, clear=False):
        if clear:
            geo.clear()

        geo.update({
            self.domain: {self.r: {"min": 0, "max": c.R}}
        })
    
    def process_parameters(self, all_params: dict, particle_params: dict, clear=False):
        if clear:
            all_params.clear()

        all_params.update(
            {key.name : value for key, value in particle_params.items()}
        )



current_param = pybamm.Parameter("Input Current") 

model = pybamm.BaseModel()
positive = SingleParticle("Positive Particle", +1)
negative = SingleParticle("Negative Particle", -1)

positive.process_model(model, current=current_param)
negative.process_model(model, current=current_param)

# def add_voltage(model: pybamm.BaseModel, positive, negative):

    # RTF = 2 * c.R_GAS * c.T / c.F
    # volmer_p = RTF * pybamm.arcsinh(positive.j / (2 * positive.j0)) 
    # volmer_n = RTF * pybamm.arcsinh(negative.j / (2 * negative.j0)) 

    # up = func(positive.surf_conc / positive. c.T)
    # v = volmer_p + up - volmer_n - un

    # """ REFERENCE
    # RT_F = param.R * T / param.F
    # j0_n = param.n.prim.j0(param.c_e_init_av, c_s_surf_n, T)
    # j0_p = param.p.prim.j0(param.c_e_init_av, c_s_surf_p, T)
    # eta_n = (2 / param.n.prim.ne) * RT_F * pybamm.arcsinh(j_n / (2 * j0_n))
    # eta_p = (2 / param.p.prim.ne) * RT_F * pybamm.arcsinh(j_p / (2 * j0_p))
    # phi_s_n = 0

    # phi_e = -eta_n - param.n.prim.U(sto_surf_n, T)
    # phi_s_p = eta_p + phi_e + param.p.prim.U(sto_surf_p, T)
    # V = phi_s_p
    # """
    # pass

# add_voltage(model, positive, negative)

geo = {}
positive.process_geometry(geo)
negative.process_geometry(geo)

param_dict = {}
positive.process_parameters(param_dict, {
    current_param: "[input]", 
    positive.conc_0: c.POS_CSN_INITIAL,
    positive.L: c.POS_ELEC_THICKNESS,
    positive.eps_n: c.POS_ELEC_POROSITY
})

negative.process_parameters(param_dict, {
    current_param: "[input]", 
    negative.conc_0: c.NEG_CSN_INITIAL,
    negative.L: c.NEG_ELEC_THICKNESS,
    negative.eps_n: c.NEG_ELEC_POROSITY
})

param = pybamm.ParameterValues(param_dict)
param.process_model(model)
param.process_geometry(geo)

PTS = 20
mesh = pybamm.Mesh(geo, 
    { positive.domain: pybamm.Uniform1DSubMesh, negative.domain: pybamm.Uniform1DSubMesh }, 
    { positive.r: PTS, negative.r: PTS }
)

disc = pybamm.Discretisation(mesh, 
    { positive.domain: pybamm.FiniteVolume(), negative.domain: pybamm.FiniteVolume() }
)
disc.process_model(model)


### SETUP DONE ###
pos_capacity = (c.POS_CSN_MAX - c.POS_CSN_INITIAL) * c.POS_ELEC_THICKNESS * (1-c.POS_ELEC_POROSITY) 
neg_capacity = (c.NEG_CSN_INITIAL - c.NEG_CSN_MIN) * c.NEG_ELEC_THICKNESS * (1-c.NEG_ELEC_POROSITY)

capacity = min(pos_capacity, neg_capacity) * (c.F / 3600) # conversion into Ah

solver = pybamm.ScipySolver()

seconds = c.RUNTIME_HOURS * 3600
time_steps = np.linspace(0, seconds, int(seconds) // 30)
print(f"Evaluating @ {len(time_steps)} timesteps")

calc_current = (capacity / c.RUNTIME_HOURS)

print(f"Discharging @ {calc_current:.3f} A, maxing electrode in {seconds} seconds")

# Evaluate concentration @ each <time_steps> steps @ at <PTS> locations from r=0->R
solution = solver.solve(model, time_steps, inputs={current_param.name: calc_current})
# solution.plot([positive.conc_name, negative.conc_name])


import pickle

fle = open("Up_func.pkl", 'rb')
Up = pickle.load(fle)
fle.close()

fle = open("Un_func.pkl", 'rb')
Un = pickle.load(fle)
fle.close()

def post_process_voltage(solution):
    global Up
    global Un

    ### NAIVE METHOD OF VOLTAGE CALCULATION (post processing)
    # Voltage = UP + VOL_P - UN - VOL_N

    voltages = []
    RTF = c.R_GAS * c.T / c.F

    surf_p = solution[positive.surf_conc_name].entries # length 24 (one at each time)
    surf_n = solution[negative.surf_conc_name].entries

    ## j (electrode current density is constant throughout?)
    j_p = solution[positive.j_name].entries[0]
    j_n = solution[negative.j_name].entries[0]

    time_steps = len(solution.t)

    for i in range(time_steps):
        inst_surf_p = surf_p[i]
        scaled_surf_p = inst_surf_p / c.POS_CSN_MAX
        j0_p = sqrt(scaled_surf_p) * sqrt(1 - scaled_surf_p)

        volmer_p = 2 * RTF * pybamm.arcsinh(j_p / (2 * j0_p))
        up = Up(scaled_surf_p)
        
        inst_surf_n = surf_n[i]
        scaled_surf_n = inst_surf_n / c.NEG_CSN_MAX
        j0_n = sqrt(scaled_surf_n) * sqrt(1 - scaled_surf_n)
        
        volmer_n = 2 * RTF * pybamm.arcsinh(j_n / (2 * j0_n))
        un = Un(scaled_surf_n)
        
        v = up + volmer_p - volmer_n - un
        voltages.append(v.value)

    return voltages

voltages = post_process_voltage(solution)

plt.plot(list(solution.t), voltages)
plt.show()

# for v in voltages:
    # print(v)