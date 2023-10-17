# at surface of particle, Neumann (we know the charge rate --> ion travel). 
# 'Molar Flux' function of current (faraday's constant converts current -> mol)
# IF concentration was known instead, THEN Dirichlet

# what's the concentration of Li @ t=0
# depends on real vs scaled concentration (latter is normalized to 0-1)

# Graphite diffusion coefficient depends on lithium concentration in that electrode
# lets assume its constant now (for modeling purpose)

import pybamm
import numpy as np
import sys

model = pybamm.BaseModel()

D = pybamm.Parameter("Diffusion Coefficient")
R = pybamm.Parameter("Particle Max Radius")
j = pybamm.Parameter("Current Density")
c0 = pybamm.Parameter("Initial Concentration")

## Faraday constant (so need not be parameter, right?) 
F = pybamm.Scalar(96485, "Faraday Constant")

c = pybamm.Variable("Concentration", domain="electrode")

# gradient is derivative in spatial dimensions
molar_flux = D * -pybamm.grad(c)    

## why is divergence tagged with negative. Makes sense for gradient due to definition of flux
dcdt = -pybamm.div(molar_flux)

model.rhs = {c: dcdt}

model.initial_conditions = {
    c: c0
}

# inner boundary condition at r=0 (center of the particle)
# no motivation for li+ concentration change -- because concentration on either side of center is same  
# flux should be 0 at that 'boundary'. So, we still use Neumann, but a constant value of 0

model.boundary_conditions = {
    c: {
        "left": (0, "Neumann"),
        "right": (-j / (F*D), "Neumann") 
    }
}

model.variables = {
    "Concentration" : c, "Flux": molar_flux
}

r = pybamm.SpatialVariable("Radius", domain="electrode", coord_sys="spherical polar")

geo = {"electrode": {r: {"min": 0, "max": R}}}

param = pybamm.ParameterValues({
    "Diffusion Coefficient": 3.9e-14,
    "Particle Max Radius": 10e-5,
    "Current Density": "[input]", #1.4
    "Initial Concentration": 2.5e4
})

param.process_model(model)
param.process_geometry(geo)

PTS = 10
mesh = pybamm.Mesh(geo, {"electrode": pybamm.Uniform1DSubMesh}, {r: PTS})
disc = pybamm.Discretisation(mesh, {"electrode": pybamm.FiniteVolume()})
disc.process_model(model)

solver = pybamm.ScipySolver()
seconds = int(sys.argv[1])
time_steps = np.linspace(0, seconds, (seconds // int(sys.argv[2])) or 100)

# Evaluate concentration @ each 100 time steps @ at 50 locations from r=0->R
solution = solver.solve(model, time_steps, inputs={"Current Density": 1.4})
solution.plot(list(model.variables.keys()))

# concs = solution['Concentration']
## returns PTS+2 entries (is this to simulate boundaries)
# print(concs(t=10))

# solution = solver.solve(model, time_steps, inputs={"Current Density": 2.8})
# solution.plot(list(model.variables.keys()))
