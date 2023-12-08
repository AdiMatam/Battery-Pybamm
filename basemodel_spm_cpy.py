import pybamm
import numpy as np
from consts import RUNTIME_HOURS

def voltage_compare(iapp: float):
    model = pybamm.lithium_ion.BasicSPM(iapp=iapp)

    param = model.default_parameter_values
    # print(param)

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
    return voltages

    # with open("voltage_compare.txt", 'w') as f:
        # voltages = solution["Voltage [V]"].entries
        # for v in voltages:
            # f.write(str(v) + '\n')

# with open("negative_concentrations.txt", 'w') as f:
    # voltages = solution["Negative particle surface concentration [mol.m-3]"].entries[0]
    # for v in voltages:
        # f.write(str(v) + '\n')

# with open("positive_concentrations.txt", 'w') as f:
    # voltages = solution["Positive particle surface concentration [mol.m-3]"].entries[0]
    # for v in voltages:
        # f.write(str(v) + '\n')

