import pybamm
import numpy as np
from matplotlib import pyplot as plt
import consts as c

from single_particle import SingleParticle

iapp = pybamm.Parameter("Input Current / Area") 
cell_voltage = "Voltage"

model = pybamm.BaseModel()
positive = SingleParticle("Positive Particle", +1, iapp)
negative = SingleParticle("Negative Particle", -1, iapp)

positive.process_model(model)
negative.process_model(model)
model.variables.update({
    cell_voltage: positive.voltage - negative.voltage
})

geo = {}
positive.process_geometry(geo)
negative.process_geometry(geo)

## TODO: improve how this is done. Parameters in a file? json or something may simplify? 
param_dict = {
    iapp.name: "[input]"
}
positive.process_parameters(param_dict, {
    positive.conc_0:    c.POS_CSN_INITIAL,
    positive.L:         c.POS_ELEC_THICKNESS,
    positive.eps_n:     c.POS_ELEC_POROSITY,
    positive.conc_max:  c.POS_CSN_MAX,
    positive.j0:        c.POS_EXCHANGE_CURRENT_DENSITY,
    positive.ocp:       c.POS_OPEN_CIRCUIT_POTENTIAL
})

negative.process_parameters(param_dict, {
    negative.conc_0:    c.NEG_CSN_INITIAL,
    negative.L:         c.NEG_ELEC_THICKNESS,
    negative.eps_n:     c.NEG_ELEC_POROSITY,
    negative.conc_max:  c.NEG_CSN_MAX,
    negative.j0:        c.NEG_EXCHANGE_CURRENT_DENSITY,
    negative.ocp:       c.NEG_OPEN_CIRCUIT_POTENTIAL
    
})

param = pybamm.ParameterValues(param_dict)
param.process_model(model)
param.process_geometry(geo)

PTS = 30
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

capacity = min(pos_capacity, neg_capacity) 
if (capacity == pos_capacity):
    print("Pos electrode parameters LIMITING")
else:
    print("Neg electrode parameters LIMITING")

capacity *= (c.F / 3600) # conversion into Ah
calc_current = (capacity / c.RUNTIME_HOURS)

seconds = c.RUNTIME_HOURS * 3600
time_steps = np.linspace(0, seconds, 250)
print(f"Evaluating @ {len(time_steps)} timesteps")
print(f"Discharging @ {calc_current:.3f} A/m2; Runtime: {seconds} seconds")

# Evaluate concentration @ each <time_steps> steps @ at <PTS> locations from r=0->R
# use NEGATIVE CURRENT <-> REPRESENTING DISCHARGE
solver = pybamm.ScipySolver()
solution = solver.solve(model, time_steps, inputs={iapp.name: -calc_current})
# solution.plot([
    # *model.variables.keys(),
    # "Voltage"
# ])

### COMPARE EXERCISE
from basemodel_spm_cpy import model_compare
my_voltages = solution[cell_voltage].entries
my_pos = solution[positive.surf_conc_name].entries[0]
my_neg = solution[negative.surf_conc_name].entries[0]
pyb_pos, pyb_neg, pyb_voltages = model_compare(calc_current)

# print(len(my_voltages))
# print(len(pyb_voltages))
# assert(len(my_voltages) == len(pyb_voltages))

plt.plot(solution.t, my_pos, label="My Model", color='r')
plt.plot(solution.t, pyb_pos, label="Pybamm Model", color='b')

plt.plot(solution.t, my_neg, label="My Model", color='r')
plt.plot(solution.t, pyb_neg, label="Pybamm Model", color='b')
plt.legend()
plt.show()

