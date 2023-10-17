import pybamm
import numpy as np
from matplotlib import pyplot as plt
import consts as c

from single_particle import SingleParticle

current_param = pybamm.Parameter("Input Current / Area") 

model = pybamm.BaseModel()
positive = SingleParticle("Positive Particle", +1, current_param)
negative = SingleParticle("Negative Particle", -1, current_param)

positive.process_model(model)
negative.process_model(model)

geo = {}
positive.process_geometry(geo)
negative.process_geometry(geo)

## TODO: improve how this is done. Parameters in a file? json or something may simplify? 
param_dict = {}
positive.process_parameters(param_dict, {
    current_param: "[input]", 
    positive.conc_0: c.POS_CSN_INITIAL,
    positive.L: c.POS_ELEC_THICKNESS,
    positive.eps_n: c.POS_ELEC_POROSITY,
    positive.conc_max: c.POS_CSN_MAX
})

negative.process_parameters(param_dict, {
    current_param: "[input]", # send value at solve-time
    negative.conc_0: c.NEG_CSN_INITIAL,
    negative.L: c.NEG_ELEC_THICKNESS,
    negative.eps_n: c.NEG_ELEC_POROSITY,
    negative.conc_max: c.NEG_CSN_MAX
})

param = pybamm.ParameterValues(param_dict)
param.process_model(model)
param.process_geometry(geo)

PTS = 50
mesh = pybamm.Mesh(geo, 
    { positive.domain: pybamm.Uniform1DSubMesh, negative.domain: pybamm.Uniform1DSubMesh }, 
    { positive.r: PTS, negative.r: PTS }
)

disc = pybamm.Discretisation(mesh, 
    { positive.domain: pybamm.FiniteVolume(), negative.domain: pybamm.FiniteVolume() }
)
disc.process_model(model)


### SETUP DONE ###
pos_capacity = (c.POS_CSN_MAX - c.POS_CSN_INITIAL) * c.POS_ELEC_THICKNESS * (1-c.POS_ELEC_POROSITY) 
neg_capacity = (c.NEG_CSN_INITIAL - c.NEG_CSN_MIN) * c.NEG_ELEC_THICKNESS * (1-c.NEG_ELEC_POROSITY)

capacity = min(pos_capacity, neg_capacity) * (c.F / 3600) # conversion into Ah

solver = pybamm.ScipySolver()

seconds = c.RUNTIME_HOURS * 3600
time_steps = np.linspace(0, seconds, int(seconds) // 60)
print(f"Evaluating @ {len(time_steps)} timesteps")

calc_current = (capacity / c.RUNTIME_HOURS)

print(f"Discharging @ {calc_current:.3f} A/m2; Runtime: {seconds} seconds")

# Evaluate concentration @ each <time_steps> steps @ at <PTS> locations from r=0->R
solution = solver.solve(model, time_steps, inputs={current_param.name: calc_current})
solution.plot([positive.conc_name, positive.j0_name, positive.sto_name])

from voltage_sim import post_process_voltage

voltages = post_process_voltage(solution, positive, negative)

plt.plot(list(solution.t), voltages)
plt.show()

