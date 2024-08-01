"""
-Create (variety of tests, C-rates, etc)
---Cycle 1 cell 
---Cycle pack of cells
*Without parameter variation.. we should see identical results.

-Processing output data in a plottable, friendly way

LATER LATER:
## -- Capacity Output --> should be plottable
      ## -- For charge capacity:, use i_intercalation
      ## -- i_sei is contributing to 'Capacity Loss'
"""


NUM_PARALLEL = 2
NUM_SERIES = 2
NUM_CYCLES = 3

BASE_CURRENT = 13.6319183090575
### ESTIMATED FROM 0.5C RATE

## input current (you can change to anything)
I_INPUT = BASE_CURRENT * NUM_PARALLEL

VOLTAGE_LOW_CUT = 2.0
VOLTAGE_HIGH_CUT =4.1

## Meshing and Discretization Parameters
HOURS = 2 
DISCRETE_PTS = 30
TIME_PTS = 100

# Data is outputted to this file.
DATA_OUTPUT = "data/2by2_3cycles_control"

#--------------------

import pybamm
from pack import Pack

model = pybamm.BaseModel()
geo = {}
parameters = {}

# I_INPUT / 20
pack = Pack(NUM_PARALLEL, NUM_SERIES, (VOLTAGE_LOW_CUT, VOLTAGE_HIGH_CUT), 1.3, model, geo, parameters)
pack.build(DISCRETE_PTS)

pack.cycler(I_INPUT, NUM_CYCLES, HOURS, TIME_PTS, output_path=DATA_OUTPUT+".csv")

import pickle

with open(DATA_OUTPUT+".pkl", 'wb') as f:
      pickle.dump(pack, f)
