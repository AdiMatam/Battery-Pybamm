import pybamm
import numpy as np
from consts import RUNTIME_HOURS

def model_compare(iapp: float):
    model = pybamm.lithium_ion.BasicSPM(iapp=iapp)

    param = model.default_parameter_values

    geo = model.default_geometry
    param.process_model(model)
    param.process_geometry(geo)
    mesh = pybamm.Mesh(geo, model.default_submesh_types, model.default_var_pts)
    disc = pybamm.Discretisation(mesh, model.default_spatial_methods)
    disc.process_model(model)

    solver = model.default_solver
    n = 250
    
    t_eval = np.linspace(0, 3600 * RUNTIME_HOURS, n)
    solution = solver.solve(model, t_eval)
    # print(solution["Negative particle surface concentration [mol.m-3]"].entries)
    # solution.plot(list(model.variables.keys()))
    # print(list(model.variables.keys()))

    voltages = solution["Voltage [V]"].entries
    neg_concs = solution["Negative particle surface concentration [mol.m-3]"].entries[0]
    pos_concs = solution["Positive particle surface concentration [mol.m-3]"].entries[0]
    return pos_concs, neg_concs, voltages

