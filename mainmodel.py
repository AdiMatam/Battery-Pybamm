NUM_SERIES = 2
NUM_PARALLEL = 10
NUM_CYCLES = 100

BASE_CURRENT = 13.6319183090575
### ESTIMATED FROM 0.5C RATE

## input current (you can change to anything)
I_INPUT = BASE_CURRENT * NUM_PARALLEL

VOLTAGE_LOW_CUT = 2.0
VOLTAGE_HIGH_CUT =4.1

## Meshing and Discretization Parameters
HOURS = 2 
DISCRETE_PTS = 50
TIME_PTS = 100

# Data is outputted to this file.
DATA_OUTPUT = "data/2by10_100cycles_porosityX"

#--------------------

import pybamm
from pack import Pack
pybamm.set_logging_level("WARNING")

model = pybamm.BaseModel()
geo = {}
parameters = {}

# I_INPUT / 20
pack = Pack(NUM_PARALLEL, NUM_SERIES, (VOLTAGE_LOW_CUT, VOLTAGE_HIGH_CUT), I_INPUT/ 10, model, geo, parameters)
pack.build(DISCRETE_PTS)

pack.cycler(I_INPUT, NUM_CYCLES, HOURS, TIME_PTS, output_path=DATA_OUTPUT+".csv")

import pickle

with open(DATA_OUTPUT+".pkl", 'wb') as f:
      pickle.dump(pack, f)
