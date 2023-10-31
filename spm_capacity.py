import pybamm
import numpy as np
from matplotlib import pyplot as plt
import consts as c

from single_particle import SingleParticle
from marquis import lico2_electrolyte_exchange_current_density_Dualfoil1998 as j0p
from marquis import graphite_electrolyte_exchange_current_density_Dualfoil1998 as j0n

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
    positive.conc_max: c.POS_CSN_MAX,
    positive.j0: j0p
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


# ### SETUP DONE ###
# pos_capacity = (c.POS_CSN_MAX - c.POS_CSN_INITIAL) * c.POS_ELEC_THICKNESS * (1-c.POS_ELEC_POROSITY) 
# neg_capacity = (c.NEG_CSN_INITIAL - c.NEG_CSN_MIN) * c.NEG_ELEC_THICKNESS * (1-c.NEG_ELEC_POROSITY)

# capacity = min(pos_capacity, neg_capacity) 
# if (capacity == pos_capacity):
    # print("Pos electrode parameters LIMITING")
# else:
    # print("Neg electrode parameters LIMITING")

# capacity *= (c.F / 3600) # conversion into Ah

# solver = pybamm.ScipySolver()

# seconds = c.RUNTIME_HOURS * 3600
# time_steps = np.linspace(0, seconds, 250)
# print(f"Evaluating @ {len(time_steps)} timesteps")

# calc_current = (capacity / c.RUNTIME_HOURS)

# print(f"Discharging @ {calc_current:.3f} A/m2; Runtime: {seconds} seconds")

# # Evaluate concentration @ each <time_steps> steps @ at <PTS> locations from r=0->R
# # use NEGATIVE CURRENT <-> REPRESENTING DISCHARGE
# solution = solver.solve(model, time_steps, inputs={current_param.name: -calc_current})
# solution.plot([positive.conc_name, positive.j0_name, positive.sto_name])

# from voltage_sim import post_process_voltage

# voltages = post_process_voltage(solution, positive, negative)

# ## COMPARE WITH PYBAMM-GENERATED

# pyb_voltages = []
# with open("compare_test.txt", 'r') as f:
    # for line in f:
        # if line:
            # pyb_voltages.append(float(line))


# plt.plot(list(solution.t), voltages, label='My Model', color='r')
# plt.plot(list(solution.t), pyb_voltages, label='Pybamm Model', color='b')
# plt.xlabel("Time (s)")
# plt.ylabel("Voltage (V)")

# plt.legend()
# plt.show()



