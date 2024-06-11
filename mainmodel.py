#-----------------
"""
Change the following simulation parameters (in all caps) and run the code.
The code will output the following at each time step:
- Concentrations of Lithium at each electrode
- Voltage across each cell
- Current in each cell

The output is sent to cycle_data.csv (Can be opened in Excel)
"""


NUM_PARALLEL = 3
NUM_SERIES = 2
NUM_CYCLES = 1
BASE_CURRENT = 1.20276592916666664

## input current (you can change to anything)
I_TOTAL = BASE_CURRENT * NUM_PARALLEL

VOLTAGE_LOW_CUT = 3.5
VOLTAGE_HIGH_CUT =3.85

## Meshing and Discretization Parameters
DISCRETE_PTS = 50
TIME_PTS = 250
RUNTIME_HOURS = 20

#--------------------


import pybamm
import params as p
from pack import Pack

i_t = pybamm.Parameter("Input Current / Area") 

model = pybamm.BaseModel()
geo = {}
parameters = {i_t.name: "[input]"}

pack = Pack(NUM_PARALLEL, NUM_SERIES, model, geo, parameters, i_t)

ivp = lambda c: (p.POS_OCP(c.pos_csn_ival / c.pos_csn_maxval))
ivn = lambda c: (p.NEG_OCP(c.neg_csn_ival / c.neg_csn_maxval))

model.initial_conditions.update({
    **{ pack.iapps[i]: -BASE_CURRENT for i in range(pack.parallel) },

    **{ pack.cells[i,j].volt: ivp(pack.cells[i,j]) - ivn(pack.cells[i,j])
            for i in range(pack.series) for j in range(pack.parallel)
      }
})

model.events += [
      pybamm.Event("Min Voltage Cutoff", pack.voltsum - VOLTAGE_LOW_CUT*pack.series),
      pybamm.Event("Max Voltage Cutoff", VOLTAGE_HIGH_CUT*pack.series - pack.voltsum),
]

print(f"Theoretical estimate: {pack.capacity}")

pack.build(DISCRETE_PTS)

variables = list(model.variables.keys())
df, caps = pack.cycler(I_TOTAL, NUM_CYCLES, RUNTIME_HOURS, TIME_PTS, variables, output_path="cycle_data.csv")

print(f"Actual: {caps}")

from plotter import plot
plot(df, pack)