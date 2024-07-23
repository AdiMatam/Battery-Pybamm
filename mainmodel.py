#-----------------
"""
Change the following simulation parameters (in all caps) and run the code.
The code will output the following at each time step:
- Concentrations of Lithium at each electrode
- Voltage across each cell
- Current in each 'string' of cells

The output is sent to cycle_data.csv (Can be opened in Excel)
"""


NUM_PARALLEL = 2
NUM_SERIES = 1
NUM_CYCLES = 50

BASE_CURRENT = 13.6319183090575 #2.4
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
# TODO: Need plotting interface so it is easier to create graphs!
DATA_OUTPUT = "capdata.csv"

#--------------------

import pybamm
from pack import Pack

model = pybamm.BaseModel()
geo = {}
parameters = {}

pack = Pack(NUM_PARALLEL, NUM_SERIES, (VOLTAGE_LOW_CUT, VOLTAGE_HIGH_CUT), model, geo, parameters)

## TODO:
## -- CC + CV Charge
## -- Capacity Output --> should be plottable
      ## -- Add equation for capacity integration at the cell-level (so that we can compare capacities)
      ## -- For charge capacity:, use i_intercalation
      ## -- i_sei is contributing to 'lost' capacity (Capacity Loss)

## CHECK (ignore this)
lhs = set([a.name for a in model.rhs.keys()]) | set([a.name for a in model.algebraic.keys()])
rhs = set([a.name for a in model.initial_conditions.keys()])
print("Determination:", len(lhs) - len(rhs))
##

pack.build(DISCRETE_PTS)

df = pack.cycler(I_INPUT, NUM_CYCLES, HOURS, TIME_PTS, output_path=DATA_OUTPUT)
