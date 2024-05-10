NUM_PARALLEL = 3
NUM_CYCLES = 1
BASE_CURRENT = 1.20276592916666664
I_TOTAL = BASE_CURRENT * NUM_PARALLEL

DISCRETE_PTS = 50
TIME_PTS = 250
RUNTIME_HOURS = 20


import pybamm
import params as p
from pack import Pack

i_t = pybamm.Parameter("Input Current / Area") 

model = pybamm.BaseModel()
geo = {}
parameters = {i_t.name: "[input]"}

pack = Pack(NUM_PARALLEL, 2, model, geo, parameters, i_t)

ivp = lambda c: (p.POS_OCP(c.pos_csn_ival / c.pos_csn_maxval))
ivn = lambda c: (p.NEG_OCP(c.neg_csn_ival / c.neg_csn_maxval))

model.initial_conditions.update({
    **{ pack.iapps[i]: -BASE_CURRENT for i in range(pack.parallel) },

    **{ pack.cells[i,j].volt: ivp(pack.cells[i,j]) - ivn(pack.cells[i,j])
            for i in range(pack.series) for j in range(pack.parallel)
      }
})

pack.build(DISCRETE_PTS)

# print(list(map(lambda x: x.name, pack.iapps)))
variables = list(model.variables.keys())

df = pack.cycler(I_TOTAL, 2, RUNTIME_HOURS, TIME_PTS, variables, output_path="cycle_data.csv")

from plotter import plot
plot(df, pack)
