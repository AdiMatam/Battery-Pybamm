import pybamm
import numpy as np
import consts as c

from single_particle import SingleParticle
from cell import Cell

"""
https://docs.pybamm.org/en/latest/source/examples/notebooks/solvers/dae-solver.html#

cell1 = Cell("Num1", iapp)
cell1.process_model(model)
cell1.process_geometry(geometry)
cell1.process_parameters(params)

cell1.voltage <- V = Ep - En

cell2 = Cell("Num2", iapp)
cell2.process_model(model)
cell2.process_geometry(geometry)
cell2.process_parameters(params)

cell2.voltage

net_voltage = pybamm.Variable("Enforced Voltage") ## this shouldn't have a domain? 

model.algebraic = {
    net_voltage: cell1.voltage - cell2.voltage
}
"""

current_param = pybamm.Parameter("Input Current / Area") 

model = pybamm.BaseModel()
geo = {}
param_dict = {}

cell1 = Cell("Num1", model, geo, current_param)
cell1.process_parameters(param_dict) ## hardcoded internally now

cell2 = Cell("Num2", model, geo, current_param)
cell2.process_parameters(param_dict) 

param = pybamm.ParameterValues(param_dict)
param.process_model(model)
param.process_geometry(geo)

# PTS = 30
# mesh = pybamm.Mesh(geo, 
    # { positive.domain: pybamm.Uniform1DSubMesh, negative.domain: pybamm.Uniform1DSubMesh }, 
    # { positive.r: PTS, negative.r: PTS }
# )

# disc = pybamm.Discretisation(mesh, 
    # { positive.domain: pybamm.FiniteVolume(), negative.domain: pybamm.FiniteVolume() }
# )
# disc.process_model(model)


# # Evaluate concentration @ each <time_steps> steps @ at <PTS> locations from r=0->R
# # use NEGATIVE CURRENT <-> REPRESENTING DISCHARGE
# calc_current = -24
# solver = pybamm.ScipySolver()
# time_steps = np.linspace(0, 3600, 250)
# solution = solver.solve(model, time_steps, inputs={current_param.name: -calc_current})