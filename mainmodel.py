### ESTIMATED FROM THEORETICAL CAPACITY FORMULATION
THEORETICAL_CAPACITY = 27.263836618115

### CHANGE SIMULATION PROFILE / OPERATING CONDITIONS HERE
# ------------------

NUM_SERIES = 5
NUM_PARALLEL = 5
NUM_CYCLES = 100

C_RATE = 1.0
BASE_CURRENT = THEORETICAL_CAPACITY * C_RATE

## input current (you can change to anything)
I_INPUT = BASE_CURRENT * NUM_PARALLEL
CURRENT_CUT_FACTOR = 1/10

VOLTAGE_LOW_CUT = 2.5 * NUM_SERIES
VOLTAGE_HIGH_CUT =4.1 * NUM_SERIES

## Meshing and Discretization Parameters
DISCRETE_PTS = 30
HOURS = 2 
TIME_PTS = 100

# Data is outputted to this file.
EXPERIMENT = "5by5_100cycles_const"

#--------------------




### DON'T CHANGE BELOW THIS!

import pybamm
from src.pack import Pack
pybamm.set_logging_level("WARNING")

model = pybamm.BaseModel()
geo = {}
parameters = {}

pack = Pack(
      EXPERIMENT, 
      NUM_PARALLEL, NUM_SERIES, I_INPUT, NUM_CYCLES, C_RATE,
      (VOLTAGE_LOW_CUT, VOLTAGE_HIGH_CUT), CURRENT_CUT_FACTOR, 
      model, geo, parameters
)

pack.export_profile()

pack.build(DISCRETE_PTS)
pack.cycler(HOURS, TIME_PTS)

