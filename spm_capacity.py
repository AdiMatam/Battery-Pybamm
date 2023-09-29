import pybamm
import numpy as np

import consts as c

class SingleParticle:

    def __init__(self, name: str, charge: int):
        self.name = name
        self.charge = charge
        self.domain = name + " dDomain"

        self.conc_0 = pybamm.Parameter(name + " pInitial Concentration")

        self.L = pybamm.Parameter(name + " pElectrode Thickness")
        self.eps_n = pybamm.Parameter(name + " pElectrode Porosity")

        self.conc_variable_name = name + " vConcentration"
        self.flux_variable_name = name + " vFlux"

        self.conc_variable = pybamm.Variable(self.conc_variable_name, domain=self.domain)
        self.r = pybamm.SpatialVariable(name + " svRadius", domain=self.domain, coord_sys="spherical polar")

    def process_model(self, model: pybamm.BaseModel, current, clear=False):
        if clear:
            model = pybamm.BaseModel()

        flux = c.D * -pybamm.grad(self.conc_variable)
        dcdt = -pybamm.div(flux)

        self.surf_conc = pybamm.surf(self.conc_variable)

        model.rhs.update({
            self.conc_variable: dcdt,
        })

        model.initial_conditions.update({
            self.conc_variable: self.conc_0,
        }) 
        
        a_term = (3 * (1 - self.eps_n)) / c.R

        self.j0 = pybamm.sqrt(self.surf_conc) * pybamm.sqrt(1 - self.surf_conc)
        self.j = current / (self.L * a_term)

        model.boundary_conditions.update({
            self.conc_variable: {
                "left":  (0, "Neumann"),
                "right": (self.j / (self.charge * c.F * c.D), "Neumann")
            },
        })
        model.variables.update({
            self.conc_variable_name: self.conc_variable, self.flux_variable_name: flux,
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

def add_voltage(model: pybamm.BaseModel, positive, negative):

    RTF = 2 * c.R_GAS * c.T / c.F
    volmer_p = RTF * pybamm.arcsinh(positive.j / (2 * positive.j0)) 
    volmer_n = RTF * pybamm.arcsinh(negative.j / (2 * negative.j0)) 

    up = func(positive.surf_conc / positive. c.T)
    v = volmer_p + up - volmer_n - un

    """
    RT_F = param.R * T / param.F
    j0_n = param.n.prim.j0(param.c_e_init_av, c_s_surf_n, T)
    j0_p = param.p.prim.j0(param.c_e_init_av, c_s_surf_p, T)
    eta_n = (2 / param.n.prim.ne) * RT_F * pybamm.arcsinh(j_n / (2 * j0_n))
    eta_p = (2 / param.p.prim.ne) * RT_F * pybamm.arcsinh(j_p / (2 * j0_p))
    phi_s_n = 0

    phi_e = -eta_n - param.n.prim.U(sto_surf_n, T)
    phi_s_p = eta_p + phi_e + param.p.prim.U(sto_surf_p, T)
    V = phi_s_p
    """
    pass

add_voltage(model, positive, negative)

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

calc_current = (capacity / c.RUNTIME_HOURS)

print(f"Discharging @ {calc_current:.3f} A, maxing electrode in {seconds} seconds")

# Evaluate concentration @ each <time_steps> steps @ at <PTS> locations from r=0->R
solution = solver.solve(model, time_steps, inputs={current_param.name: calc_current})
print(solution.t)
print(solution[positive.conc_variable_name].entries)
print(solution[negative.conc_variable_name].entries)

solution.plot(list(model.variables.keys()))