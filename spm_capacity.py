import pybamm
import numpy as np
import sys

# assuming particles behave identically. Only boundary conditions are different

I = pybamm.Parameter("Input Current / Area") 

# assuming "Diffusion Coefficient" is constant w/ respect to Concentration @ r
D = 3.9e-14 # pybamm.Parameter("pDiffusion Coefficient")
F = 96485

# not scaling radius at the moment... 
R = 5.5e-06 

class SingleParticle:

    def __init__(self, name: str, charge: int):
        self.name = name
        self.charge = charge
        self.domain = name + " dDomain"

        self.conc_0 = pybamm.Parameter(name + " pInitial Concentration")

        self.L = pybamm.Parameter(name + "pElectrode Thickness")
        self.eps_n = pybamm.Parameter(name + "pElectrode Porosity")

        self.conc = pybamm.Variable(name + " vConcentration", domain=self.domain)
        self.r = pybamm.SpatialVariable(name + " svRadius", domain=self.domain, coord_sys="spherical polar")

    def process_model(self, model: pybamm.BaseModel, clear=False):
        if clear:
            model = pybamm.BaseModel()


        flux = D * -pybamm.grad(self.conc)
        dcdt = -pybamm.div(flux)

        model.rhs.update({
            self.conc: dcdt
        })

        model.initial_conditions.update({
            self.conc: self.conc_0
        }) 
        
        model.boundary_conditions.update({
            self.conc: {
                "left":  (0, "Neumann"),
                "right": (I / (self.charge * F * D * self.L * ((3 * (1 - self.eps_n)) / R)), "Neumann")
            },
        })
        model.variables.update({
            self.conc.name: self.conc, f"{self.name} Flux": flux
        })
    
    def process_geometry(self, geo: dict, clear=False):
        if clear:
            geo.clear()

        geo.update({
            self.domain: {self.r: {"min": 0, "max": R}}
        })
    
    def process_parameters(self, params: dict, model, geo: dict):
        create_params = {key.name : value for key, value in params.items()}
        param = pybamm.ParameterValues(create_params)
        param.process_model(model)
        param.process_geometry(geometry)


model = pybamm.BaseModel()
positive = SingleParticle("Positive Particle", +1)

positive.process_model(model)

geometry = {}
positive.process_geometry(geometry)

positive.process_parameters({
    I: 5,
    positive.conc_0: 2.5e4,
    positive.L: 7.56e-05,
    positive.eps_n: 0.33
}, model, geometry)


PTS = 10
mesh = pybamm.Mesh(geometry, 
    {positive.domain: pybamm.Uniform1DSubMesh}, 
    {positive.r: PTS}
)

disc = pybamm.Discretisation(mesh, 
    {positive.domain: pybamm.FiniteVolume()}
)
disc.process_model(model)

# solve
solver = pybamm.ScipySolver()
seconds = int(sys.argv[1])
time_steps = np.linspace(0, seconds, (seconds // int(sys.argv[2])) or 100)

# Evaluate concentration @ each <time_steps> steps @ at 10 locations from r=0->1
solution = solver.solve(model, time_steps)
solution.plot(list(model.variables.keys()))