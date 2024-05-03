NUM_PARALLEL = 3
NUM_CYCLES = 1
BASE_CURRENT = 1.20276592916666664
I_TOTAL = BASE_CURRENT * NUM_PARALLEL
VOLTAGE_CUTOFF = (2.0, 5.0) # effectively disabled

DISCRETE_PTS = 50
TIME_PTS = 250
RUNTIME_HOURS = 20


import pybamm
import numpy as np
from consts import F, R_GAS, T
import params as p
from pack import Pack
from cell import Cell

i_t = pybamm.Parameter("Input Current / Area") 

model = pybamm.BaseModel()
geo = {}
parameters = {i_t.name: -I_TOTAL}

pack = Pack(NUM_PARALLEL, 2, model, geo, parameters, i_t, VOLTAGE_CUTOFF)

ivp = lambda c: (p.POS_OCP(c.pos_csn_ival / c.pos_csn_maxval))
ivn = lambda c: (p.NEG_OCP(c.neg_csn_ival / c.neg_csn_maxval))

model.initial_conditions.update({
    **{ pack.iapps[i]: -BASE_CURRENT for i in range(pack.parallel) },

    **{ pack.cells[i,j].volt: ivp(pack.cells[i,j]) - ivn(pack.cells[i,j])
            for i in range(pack.series) for j in range(pack.parallel)
      }
})

for value in model.initial_conditions.values():
    print(value)

pack.build(DISCRETE_PTS)

solver = pybamm.CasadiSolver()
time_steps = np.linspace(0, 3600 * RUNTIME_HOURS, TIME_PTS)

solution = solver.solve(model, time_steps)
for t in time_steps:
    # print(solution[c1.pos.surf_csn_name](t)[-1])
    print(solution[pack.cells[0,0].volt.name](t) + solution[pack.cells[1,0].volt.name](t))
    print(solution[pack.cells[0,1].volt.name](t) + solution[pack.cells[1,1].volt.name](t))
    print()

quit()