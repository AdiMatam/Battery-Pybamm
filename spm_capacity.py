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

        self.conc = pybamm.Variable(name + " vConcentration", domain=self.domain)
        self.r = pybamm.SpatialVariable(name + " svRadius", domain=self.domain, coord_sys="spherical polar")

    def process_model(self, model: pybamm.BaseModel, current, clear=False):
        if clear:
            model = pybamm.BaseModel()

        flux = c.D * -pybamm.grad(self.conc)
        dcdt = -pybamm.div(flux)

        model.rhs.update({
            self.conc: dcdt
        })

        model.initial_conditions.update({
            self.conc: self.conc_0
        }) 
        
        a_term = (3 * (1 - self.eps_n)) / c.R
        model.boundary_conditions.update({
            self.conc: {
                "left":  (0, "Neumann"),
                "right": (current / (self.charge * c.F * c.D * self.L * a_term), "Neumann")
            },
        })
        model.variables.update({
            self.conc.name: self.conc, f"{self.name} Flux": flux
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
time_steps = np.linspace(0, seconds, 60)

# Evaluate concentration @ each <time_steps> steps @ at <PTS> locations from r=0->R
calc_current = (capacity / c.RUNTIME_HOURS)

print(f"Discharging @ {calc_current:.3f} A, maxing electrode in {seconds} seconds")

solution = solver.solve(model, time_steps, inputs={current_param.name: calc_current})
solution.plot(list(model.variables.keys()))