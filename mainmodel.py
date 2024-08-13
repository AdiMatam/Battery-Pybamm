"""
A given simulation run needs a 'operating condition' attribute
--> Make sure Temperature, time, C-rate all stored
--> Param values are stored, but NOT config. Make sure latter is also stored!! 
----> VIA Profile object ??
"""


NUM_SERIES = 1
NUM_PARALLEL = 2
NUM_CYCLES = 5

# 13.6319183090575
BASE_CURRENT = 13.3
### ESTIMATED FROM 0.5C RATE

## input current (you can change to anything)
I_INPUT = BASE_CURRENT * NUM_PARALLEL

VOLTAGE_LOW_CUT = 2.5
VOLTAGE_HIGH_CUT =4.1

## Meshing and Discretization Parameters
HOURS = 2 
DISCRETE_PTS = 10
TIME_PTS = 100

# Data is outputted to this file.
EXPERIMENT = "1by2_5cycles_isocz"

#--------------------

import pybamm
from pack import Pack
pybamm.set_logging_level("WARNING")

model = pybamm.BaseModel()
geo = {}
parameters = {}

pack = Pack(EXPERIMENT, NUM_PARALLEL, NUM_SERIES, I_INPUT, NUM_CYCLES, (VOLTAGE_LOW_CUT, VOLTAGE_HIGH_CUT), 1/10, model, geo, parameters)

pack.export_profile()

# pack.build(DISCRETE_PTS)

# pack.cycler(HOURS, TIME_PTS)

# import pickle
# with open(f"data/{EXPERIMENT}/model.pkl", 'wb') as f:
      # pickle.dump(pack, f)
