"""
Heat Equation: Tt = kTxx + Q(x)

u(0, t) = 0
u(L, t) = 0
u(x, 0) = 2x - x^2

Q(x) = 1 - |x - 1|
"""

import pybamm
import numpy as np

model = pybamm.BaseModel()

x = pybamm.SpatialVariable("x", domain="rod", coord_sys="cartesian")
u = pybamm.Variable("Temperature", domain="rod")
k = pybamm.Parameter("Thermal diffusivity")

# kTxx
# GRADIENT: in scalar-valued -> out vector-valued
# partial derivative with respect to each spatial direction (only 'x' in this case)
# vector-valued function --> generates a 'vector field' indicating partial du/dx at each x
flux = -k * pybamm.grad(u)

# DIVERGENCE: in vector-valued -> out scalar-valued
### WHY DO WE USE DIVERGENCE? SHOULDN'T WE DO ANOTHER GRADIENT TO FIND SECOND DERIV?
### LOOK AT KHAN ACADEMY (makes sense)
flux_deriv = -pybamm.div(flux)

# source of energy
q = 5 # 1 - pybamm.Function(np.abs, x - 1)

dudt = flux_deriv + q

# key -> variable being evaluated.
model.rhs = {u: dudt}

# at surface of particle, Neumann (we know the charge rate --> ion travel). 'Molar Flux' function of current
# IF concentration was known instead, THEN Dirichlet

# second boundary condition at r=0 (center of the particle)
# no motivation for li+ concentration change -- because concentration on either side of center is same  
# flux should be 0 at that 'boundary'. So, we still use Neumann, but a constant value of 0
model.boundary_conditions = {
    u: {
        "left": (pybamm.Scalar(0), "Dirichlet"), # temperature (or whatever value) at boundary is known 
        "right": (pybamm.Scalar(0), "Dirichlet")
    }
}
# with Neumann, given NOT value, but RATE at which thermal energy being put it into system
# heat flux is in terms of AREA (not single dimension)

# what's the concentration of Li @ t=0
# depends on real vs scaled concentration (latter is normalized to 0-1)
model.initial_conditions = {u: 2 * x - x**2}
model.variables = {'Temperature': u, 'Flux': flux, 'Source': q}


geometry = {"rod": {x: {"min": pybamm.Scalar(0), "max": pybamm.Scalar(2)}}}

# STROng function of temp / concentration
param = pybamm.ParameterValues({"Thermal diffusivity": 0.5}) #  equivalent to diffusion coeffcicent in SPM. 
param.process_model(model)
param.process_geometry(geometry)

submesh_types = {"rod": pybamm.Uniform1DSubMesh}

# number of points to 'track' over the rod geometry (so split 0< x <2 range into 30)
var_pts = {x: 30}

## slicing it into number of points
mesh = pybamm.Mesh(geometry, submesh_types, var_pts)
spatial_methods = {"rod": pybamm.FiniteVolume()}
disc = pybamm.Discretisation(mesh, spatial_methods)
disc.process_model(model)

solver = pybamm.ScipySolver()

# 100 time intervals between 0 and 1s. Simulate EACH of 
t = np.linspace(0, 1, 100)
solution = solver.solve(model, t)

# solution.t --> the linspace array
# solution.y --> Array<Array<int, TimeSteps>, Pts>
# Simulate EACH point in rod at EACH time step

print(solution['Temperature'].data)
# solution.plot(list(model.variables.keys()))
print(len(solution['Temperature'](t=0.5)))
