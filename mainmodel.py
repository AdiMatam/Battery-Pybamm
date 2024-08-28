### CHANGE SIMULATION PROFILE / OPERATING CONDITIONS HERE
# ------------------

NUM_SERIES = 5
NUM_PARALLEL = 5
NUM_CYCLES = 100

C_RATE = 1.0
# BASE_CURRENT = THEORETICAL_CAPACITY * C_RATE
# I_INPUT = BASE_CURRENT * NUM_PARALLEL

VOLTAGE_WINDOW = (
      2.5 * NUM_SERIES,
      4.1 * NUM_SERIES
)

CURRENT_CUT_FACTOR = 1/10
CAPACITY_CUT_FACTOR = 0.95

## Meshing and Discretization Parameters
DISCRETE_PTS = 30
HOURS = 2 
TIME_PTS = 100

# Data is outputted to this file.
EXPERIMENT = "5by5_100cycles_por3"

#--------------------



### DON'T CHANGE BELOW THIS!

import pybamm
from src.pack import Pack
pybamm.set_logging_level("WARNING")

model = pybamm.BaseModel()
geo = {}
parameters = {}

pack = Pack(EXPERIMENT, NUM_PARALLEL, NUM_SERIES, model, geo, parameters)
pack.set_charge_protocol(NUM_CYCLES, C_RATE, c_rate=True)
pack.set_cutoffs(VOLTAGE_WINDOW, CURRENT_CUT_FACTOR, CAPACITY_CUT_FACTOR)

pack.export_profile()

pack.build(DISCRETE_PTS)
pack.cycler(HOURS, TIME_PTS)

