from consts import THEORETICAL_CAPACITY

### CHANGE SIMULATION PROFILE / OPERATING CONDITIONS HERE
# ------------------
NUM_SERIES = 1
NUM_PARALLEL = 1
NUM_CYCLES = 30

### disable this flag and use I_INPUT to directly apply desired current
USE_C_RATE = True
C_RATE = 1.0
BASE_CURRENT = THEORETICAL_CAPACITY * C_RATE
I_INPUT = BASE_CURRENT * NUM_PARALLEL


VOLTAGE_WINDOW = (
      3.0 * NUM_SERIES,
      4.2 * NUM_SERIES
)

CURRENT_CUT_FACTOR = 1/10
CAPACITY_CUT_FACTOR = 0.80

## Meshing and Discretization Parameters
### Change 'hours' for lower/higher simulation runtime cutoff
### Change 'time_pts' for more/fewer time outputs
HOURS = (1./C_RATE) * 2.0 
TIME_PTS = 100
DISCRETE_PTS = 50

# Data is outputted to this subfolder of 'data/'.
EXPERIMENT = "Single_1.0C_3.0V"

#--------------------



### DON'T CHANGE BELOW THIS!

import pybamm
from src.pack import Pack
pybamm.set_logging_level("WARNING")

model = pybamm.BaseModel()
geo = {}
parameters = {}

pack = Pack(EXPERIMENT, NUM_PARALLEL, NUM_SERIES, model, geo, parameters)
if USE_C_RATE:
      pack.set_charge_protocol(NUM_CYCLES, C_RATE, use_c_rate=True)
else:
      pack.set_charge_protocol(NUM_CYCLES, I_INPUT, use_c_rate=False)
pack.set_cutoffs(VOLTAGE_WINDOW, CURRENT_CUT_FACTOR, CAPACITY_CUT_FACTOR)

pack.export_profile()

pack.build(DISCRETE_PTS)
pack.cycler(HOURS, TIME_PTS)

