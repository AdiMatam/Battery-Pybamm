C_RATE = 1/20 
RUNTIME_HOURS = 1 / C_RATE # hours
NUM_CELLS = 4
DISCHARGE_CURRENT = -1.2


import pybamm
import numpy as np
from cell import Cell
import params as p

model = pybamm.BaseModel()
geo = {}
i_total = pybamm.Parameter("Input Current / Area") 
#voltage = pybamm.Variable("String Voltage")

parameters = {
    i_total.name: DISCHARGE_CURRENT
}

cells = [Cell(f"Cell {i + 1}", model, geo, parameters) for i in range(NUM_CELLS)]

for i in range(len(cells) - 1):
    model.algebraic[cells[i].iapp] = cells[i + 1].iapp - (cells[i].iapp)

model.algebraic[cells[-1].iapp] = i_total - (cells[-1].iapp) # all same current

#model.algebraic[voltage] = 0
vval = 0

# for cell in cells:
    # model.algebraic[voltage] += (cell.pos.bv_term - cell.neg.bv_term)
    # vval += (p.POS_OCP(cell.pos_csn_ival / cell.pos_csn_maxval) - p.NEG_OCP(cell.neg_csn_ival / cell.neg_csn_maxval))

model.initial_conditions.update({
    #voltage: vval.value,
    **{ cell.iapp: DISCHARGE_CURRENT for cell in cells }
})

# model.variables.update({
    # voltage.name: voltage
# })

param_ob = pybamm.ParameterValues(parameters)
param_ob.process_model(model)
param_ob.process_geometry(geo)

PTS = 30
particles = [] 
for cell in cells:
    particles.append(cell.pos)
    particles.append(cell.neg)

mesh = pybamm.Mesh(geo, 
    { p.domain: pybamm.Uniform1DSubMesh for p in particles },
    { p.r: PTS for p in particles }
)

disc = pybamm.Discretisation(mesh, 
    { p.domain: pybamm.FiniteVolume() for p in particles }
)
disc.process_model(model)

solver = pybamm.CasadiSolver()
time_steps = np.linspace(0, 3600 * RUNTIME_HOURS, 250)
solution = solver.solve(model, time_steps)

