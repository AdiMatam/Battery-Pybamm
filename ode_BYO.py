"""
dx/dt = 5x + y;  x(0) = 3
dy/dt = 3x - 4y; y(0) = 25

To verify solution of coupled system by 'hand':

https://www.youtube.com/watch?v=z3Ag8WF5M_c0
https://tutorial.math.lamar.edu/classes/de/SolutionsToSystems.aspx

"""

import pybamm
import numpy as np

model = pybamm.BaseModel()
x = pybamm.Variable('x')
y = pybamm.Variable('y')

model.rhs = {
    x: 5*x + y,
    y: 3*x - 4*y
}

model.initial_conditions = {
    x: 3,
    y: 25
}

model.variables = {'x': x, 'y': y}

# No boundary conditions (we are not evaluating this equation in context of a physical system/mesh)
disc = pybamm.Discretisation()
disc.process_model(model)

solver = pybamm.ScipySolver()
t_steps = np.linspace(0, 1, 50)
solution = solver.solve(model, t_steps)

print(solution['x'].data)
print()
print(solution['y'].data)

solution.plot(['x', 'y'])
