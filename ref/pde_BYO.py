"""
Linear diffusion on unit sphere. (delta concentration dt)
This will require spatial domain.

dc/dt = dot(DEL, DEL(c)) ???

grad(SclFunction) = DEL(SclFunction)         ==> VecFunction
div (VecFunction) = dot(DEL, VecFunction)    ==> SclFunction

soo... 
dc/dt = dot(DEL, DEL(c)) ==> dot(DEL, grad(c)) ==> div(grad(c))


@ r=0
dc/dr = 0

@ r=1
dc/dr = 2

@ t=0
c = 1
"""

import pybamm
import numpy as np

model = pybamm.BaseModel()
c = pybamm.Variable('Concentration', domain="sphere")
r = pybamm.SpatialVariable("Radius", domain="sphere", coord_sys="spherical polar")
bigR = pybamm.Parameter("Max Radius")
outerBoundary = pybamm.Parameter("Outer Boundary Flux")

# definition of flux requires the negative
# if dc/dr > 0 as r -> up, then flux points the other way
dcdt = -pybamm.div(-pybamm.grad(c))

model.rhs = {c: dcdt}

model.initial_conditions = {
    c: 1
}

# Neumann conditions because we know boundary RATES not VALUES
model.boundary_conditions = {
    c: {
        "left":  (0, "Neumann"),
        "right": (outerBoundary, "Neumann")
    }
}

model.variables = {
    "Concentration": c
}

geo = {'sphere': 
    { r: {"min": 0, "max": bigR} } 
}

param = pybamm.ParameterValues({"Max Radius" : 5, "Outer Boundary Flux" : 2})
param.process_geometry(geo)
param.process_model(model)

mesh = pybamm.Mesh(
    geo, {"sphere": pybamm.Uniform1DSubMesh}, {r: 10}
)
disc = pybamm.Discretisation(mesh, {"sphere": pybamm.FiniteVolume()})

disc.process_model(model)

solver = pybamm.ScipySolver()
time_steps = np.linspace(0, 1, 100)

# Evaluate concentration @ each 100 time steps @ at 10 locations from r=0->1
solution = solver.solve(model, time_steps)

concs = solution['Concentration']
print(concs(t=0.5)) # Concentrations across radial domain @ t= 0.5
print(concs(x=0.5)) # Concentrations across time domain @ r=0.5 (not sure why parameter is 'x')
solution.plot(["Concentration"])

