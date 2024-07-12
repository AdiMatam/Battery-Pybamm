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
# BASE_CURRENT = 1.20276592916666664

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
pack = Pack(NUM_PARALLEL, NUM_SERIES, model, geo, parameters, i_input)

model.initial_conditions.update({
    **{ pack.iapps[i]: -BASE_CURRENT for i in range(pack.parallel) },
})

model.events += [
      pybamm.Event("Min Voltage Cutoff", pack.voltsum - VOLTAGE_LOW_CUT*pack.series),
      pybamm.Event("Max Voltage Cutoff", VOLTAGE_HIGH_CUT*pack.series - pack.voltsum),
]

## CHECK
lhs = set([a.name for a in model.rhs.keys()]) | set([a.name for a in model.algebraic.keys()])
rhs = set([a.name for a in model.initial_conditions.keys()])

# print(len(lhs) - len(rhs))

pack.build(DISCRETE_PTS)

df, caps = pack.cycler(I_TOTAL, NUM_CYCLES, HOURS, TIME_PTS, output_path=DATA_OUTPUT)

#print(f"Actual: {caps}")

# from ref.plotter import plot
# plot(df, pack)
