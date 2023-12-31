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
model.variables.update({
    "Voltage": positive.voltage - negative.voltage
})

model.events += [
    pybamm.Event("Voltage Min Cutoff", model.variables["Voltage"] - 3.0)
]

geo = {}
positive.process_geometry(geo)
negative.process_geometry(geo)

## TODO: improve how this is done. Parameters in a file? json or something may simplify? 
param_dict = {
    current_param.name: "[input]"
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
solution = solver.solve(model, time_steps, inputs={current_param.name: -calc_current})
solution.plot([
    *model.variables.keys(),
    "Voltage"
])