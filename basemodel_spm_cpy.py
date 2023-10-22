import pybamm
import numpy as np

model = pybamm.lithium_ion.BasicSPM()

param = model.default_parameter_values

geo = model.default_geometry
param.process_model(model)
param.process_geometry(geo)
mesh = pybamm.Mesh(geo, model.default_submesh_types, model.default_var_pts)
disc = pybamm.Discretisation(mesh, model.default_spatial_methods)
disc.process_model(model)


solver = model.default_solver
n = 250
t_eval = np.linspace(0, 3600, n)
solution = solver.solve(model, t_eval)
solution.plot(list(model.variables.keys()))

# with open("compare_test.txt", 'w') as f:
    # voltages = solution["Voltage [V]"].entries
    # for v in voltages:
        # f.write(str(v) + '\n')
