from consts import THEORETICAL_CAPACITY

### CHANGE SIMULATION PROFILE / OPERATING CONDITIONS HERE
# ------------------
NUM_SERIES = 1
NUM_PARALLEL = 1
NUM_CYCLES = 100

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
CAPACITY_CUT_FACTOR = 0.50

## Meshing and Discretization Parameters
### Change 'hours' for lower/higher simulation runtime cutoff
### Change 'time_pts' for more/fewer time outputs
HOURS = (1./C_RATE) * 2.0 
TIME_PTS = 100
DISCRETE_PTS = 100

# Data is outputted to this subfolder of 'data/'.
EXPERIMENT = "100_noage_cccv_CONTEXTSW"

#--------------------



### DON'T CHANGE BELOW THIS!

import pybamm
from pack import Pack
pybamm.set_logging_level("WARNING")

pack = Pack(EXPERIMENT, NUM_PARALLEL, NUM_SERIES)
if USE_C_RATE:
      pack.set_charge_protocol(NUM_CYCLES, C_RATE, use_c_rate=True)
else:
      pack.set_charge_protocol(NUM_CYCLES, I_INPUT, use_c_rate=False)
pack.set_cutoffs(VOLTAGE_WINDOW, CURRENT_CUT_FACTOR, CAPACITY_CUT_FACTOR)

# pack.init()
# pybamm.step.current(1, duration="1 hour", termination="2.5 V")
# experiment = pybamm.Experiment([pybamm.step.current(I_INPUT, duration=f"{HOURS} hours", termination="3.0V")])
# sim = pybamm.Simulation(pack.model, experiment=experiment)
# sim.solve()
# sim.plot()

# pack.build(DISCRETE_PTS)
pack.cycler(HOURS, TIME_PTS)

