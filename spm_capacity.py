import pybamm
import numpy as np
import consts as c

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
param_dict = {
    current_param.name: "[input]"
}

net_voltage = pybamm.Variable("Enforced Voltage") ## this shouldn't have a domain? 

cell1 = Cell("Num1", model, geo, current_param)
cell2 = Cell("Num2", model, geo, current_param)
CELLS = [cell1, cell2]

model.algebraic = {
    net_voltage: cell1.voltage - cell2.voltage
}

cell1.set_parameters(param_dict) ## hardcoded internally now
cell2.set_parameters(param_dict) 

# THE 'TOP LEVEL' CAN STILL BE ABSTRACTED? 
param = pybamm.ParameterValues(param_dict)
param.process_model(model)
param.process_geometry(geo)

PTS = 30
mesh = pybamm.Mesh(geo, 
    cell1.get_meshed_objects() | cell2.get_meshed_objects(),
    cell1.get_radial_mesh(PTS) | cell2.get_radial_mesh(PTS)
)

disc = pybamm.Discretisation(mesh, 
    cell1.get_discretized() | cell2.get_discretized()
)
disc.process_model(model)


# # Evaluate concentration @ each <time_steps> steps @ at <PTS> locations from r=0->R
# # use NEGATIVE CURRENT <-> REPRESENTING DISCHARGE
# calc_current = -24
# solver = pybamm.ScipySolver()
# time_steps = np.linspace(0, 3600, 250)
# solution = solver.solve(model, time_steps, inputs={current_param.name: -calc_current})