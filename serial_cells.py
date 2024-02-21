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
voltage = pybamm.Variable("String Voltage")

parameters = {
    i_total.name: DISCHARGE_CURRENT
}

cells = [Cell(f"Cell {i + 1}", model, geo, parameters) for i in range(NUM_CELLS)]

for i in range(len(cells) - 1):
    model.algebraic[cells[i].iapp] = cells[i + 1].iapp - (cells[i].iapp)

# model.algebraic[cells[0].pos.phi] = cells[0].pos.bv_term - (cells[0].pos.phi)
# model.algebraic[cells[-1].neg.phi] = cells[-1].neg.bv_term - (cells[-1].neg.phi)
model.algebraic[cells[-1].iapp] = i_total - (cells[-1].iapp) # all same current

model.algebraic[voltage] = 0
for cell in cells:
    model.algebraic[voltage] += (cells[0].pos.bv_term - cells[0].neg.bv_term)

pos_phi_init = p.POS_OCP(cells[0].pos_csn_ival / cells[0].pos_csn_maxval)
neg_phi_init = p.NEG_OCP(cells[-1].neg_csn_ival / cells[-1].neg_csn_maxval)

model.initial_conditions.update({
    cells[0].pos.phi: 4 * pos_phi_init,
    cells[-1].neg.phi: 4 * neg_phi_init,
    **{ cell.iapp: DISCHARGE_CURRENT for cell in cells }
})

model.variables.update({
    cells[0].pos.phi_name: cells[0].pos.phi,
    cells[-1].neg.phi_name: cells[-1].neg.phi,
})

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

volt = cells[0].pos.phi - cells[-1].neg.phi
solution.plot([cells[0].pos.phi_name, cells[-1].neg.phi_name])

### SOLVED
# from matplotlib import pyplot as plt

# fig, axs = plt.subplots(1, 4, figsize=(13,13))

# for i in range(len(axs)):
    # ax = axs[i]
    # iapp = solution[cells[i].iapp.name].entries
    # ax.plot(solution.t, np.round(iapp, 5))
    # ax.set_xlabel("t")
    # ax.set_ylabel("I (A)")
    # ax.set_title(cells[i].iapp.name)

# plt.tight_layout()
# plt.show()