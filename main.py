import pybamm

model = pybamm.lithium_ion.SPM()

param = model.default_parameter_values
print(param)

# sim = pybamm.Simulation(model)
# sim.solve([0, 3600])
# sim.plot()

# print()
