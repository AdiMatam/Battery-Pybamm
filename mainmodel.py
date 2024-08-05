NUM_SERIES = 2
NUM_PARALLEL = 5
NUM_CYCLES = 300

# 13.6319183090575
BASE_CURRENT = 13.3
### ESTIMATED FROM 0.5C RATE

## input current (you can change to anything)
I_INPUT = BASE_CURRENT * NUM_PARALLEL

VOLTAGE_LOW_CUT = 2.5
VOLTAGE_HIGH_CUT =4.1

## Meshing and Discretization Parameters
HOURS = 2 
DISCRETE_PTS = 100
TIME_PTS = 250

# Data is outputted to this file.
EXPERIMENT = "2by5_300cycles_isoc"

#--------------------

import pybamm
from pack import Pack
pybamm.set_logging_level("WARNING")

model = pybamm.BaseModel()
geo = {}
parameters = {}

pack = Pack(EXPERIMENT, NUM_PARALLEL, NUM_SERIES, (VOLTAGE_LOW_CUT, VOLTAGE_HIGH_CUT), I_INPUT/ 10, model, geo, parameters)
pack.build(DISCRETE_PTS)

pack.cycler(I_INPUT, NUM_CYCLES, HOURS, TIME_PTS)

import pickle
with open(f"data/{EXPERIMENT}/model.pkl", 'wb') as f:
      pickle.dump(pack, f)
