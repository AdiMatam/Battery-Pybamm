import pybamm
import numpy as np
import sys

import consts as c

out_current = pybamm.Parameter("Input Current") 

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

    def process_model(self, model: pybamm.BaseModel, clear=False):
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
                "right": (out_current / (self.charge * c.F * c.D * self.L * a_term), "Neumann")
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
    
    def process_parameters(self, params: dict, model, geo: dict):
        create_params = {key.name : value for key, value in params.items()}
        param = pybamm.ParameterValues(create_params)
        param.process_model(model)
        param.process_geometry(geo)


model = pybamm.BaseModel()
positive = SingleParticle("Positive Particle", +1)

positive.process_model(model)

geometry = {}
positive.process_geometry(geometry)

positive.process_parameters({
    out_current: "[input]", 
    positive.conc_0: c.POS_CSN_INITIAL,
    positive.L: c.ELEC_THICKNESS,
    positive.eps_n: c.ELEC_POROSITY
}, model, geometry)

PTS = 20
mesh = pybamm.Mesh(geometry, 
    { positive.domain: pybamm.Uniform1DSubMesh }, 
    { positive.r: PTS }
)

disc = pybamm.Discretisation(mesh, 
    { positive.domain: pybamm.FiniteVolume() }
)
disc.process_model(model)


### SETUP DONE ###
capacity = (c.POS_CSN_MAX - c.POS_CSN_INITIAL) * c.ELEC_THICKNESS * (1-c.ELEC_POROSITY) 
capacity *= (c.F / 3600) # A*h

solver = pybamm.ScipySolver()

if len(sys.argv) > 1:
    seconds = int(sys.argv[1])
    time_steps = np.linspace(0, seconds, seconds // int(sys.argv[2]))
else:
    seconds = c.RUNTIME_HOURS * 3600 # int(sys.argv[1])
    time_steps = np.linspace(0, seconds, 60)

# Evaluate concentration @ each <time_steps> steps @ at <PTS> locations from r=0->R
calc_current = (capacity / c.RUNTIME_HOURS)

print(f"Discharging @ {calc_current:.3f} A, maxing electrode in {seconds} seconds")

solution = solver.solve(model, time_steps, inputs={out_current.name: calc_current})
solution.plot(list(model.variables.keys()))
