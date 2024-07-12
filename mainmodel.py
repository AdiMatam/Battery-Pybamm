#-----------------
"""
Change the following simulation parameters (in all caps) and run the code.
The code will output the following at each time step:
- Concentrations of Lithium at each electrode
- Voltage across each cell
- Current in each cell

The output is sent to cycle_data.csv (Can be opened in Excel)
"""


NUM_PARALLEL = 2
NUM_SERIES = 1
NUM_CYCLES = 1

BASE_CURRENT = 13.6319183090575 #2.4
## input current (you can change to anything)
I_TOTAL = BASE_CURRENT * NUM_PARALLEL

VOLTAGE_LOW_CUT = 2.0
VOLTAGE_HIGH_CUT =4.1

## Meshing and Discretization Parameters
HOURS = 2 
DISCRETE_PTS = 30
TIME_PTS = 100

# Data is outputted to this file.
# TODO: I will create a plotting interface so it is easier to plot different characteristics!
DATA_OUTPUT = "mydata.csv"

#--------------------

import pybamm
from pack import Pack

i_input = pybamm.Parameter("Input Current / Area") 

model = pybamm.BaseModel()
geo = {}
parameters = {i_input.name: "[input]"}

# TODO: Iapp should be first argument.
pack = Pack(i_input, NUM_PARALLEL, NUM_SERIES, (VOLTAGE_LOW_CUT, VOLTAGE_HIGH_CUT), model, geo, parameters)

## CHECK
lhs = set([a.name for a in model.rhs.keys()]) | set([a.name for a in model.algebraic.keys()])
rhs = set([a.name for a in model.initial_conditions.keys()])

print("Determination:", len(lhs) - len(rhs))

pack.build(DISCRETE_PTS)

df, caps = pack.cycler(I_TOTAL, NUM_CYCLES, HOURS, TIME_PTS, output_path=DATA_OUTPUT)

import pickle

with open("cells.pkl", "wb") as f:
      pickle.dump(pack.cells, f)

