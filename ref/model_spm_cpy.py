import pybamm
import numpy as np
import pickle


model = pybamm.lithium_ion.SPM()

# print(model.variables.keys())
param = model.default_parameter_values
print(param)

geo = model.default_geometry
param.process_model(model)
param.process_geometry(geo)
mesh = pybamm.Mesh(geo, model.default_submesh_types, model.default_var_pts)
disc = pybamm.Discretisation(mesh, model.default_spatial_methods)
disc.process_model(model)


solver = model.default_solver
n = 250
t_eval = np.linspace(0, 3600, n)
# print('Solving using',type(solver).__name__,'solver...')
solution = solver.solve(model, t_eval)
# print('Finished.')

a = solution['Positive electrode open-circuit potential [V]'].entries
#print(a)
# print(len(a[0]))

# param = model.default_parameter_values
# print(param)
pos_ocp = param['Positive electrode OCP [V]']
neg_ocp = param['Negative electrode OCP [V]']

# fle = open("Un_func.pkl", 'wb')
# pickle.dump(neg_ocp, fle)
# fle.close()

# fle = open("Up_func.pkl", 'wb')
# pickle.dump(pos_ocp, fle)
# fle.close()

print(neg_ocp)
print(type(pos_ocp))

# sim = pybamm.Simulation(model)
# sim.solve([0, 3600])
# sim.plot()

# print()
