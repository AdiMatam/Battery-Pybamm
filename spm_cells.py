C_RATE = 1/20 
RUNTIME_HOURS = 1 / C_RATE # hours
NUM_CELLS = 4
DISCHARGE_CURRENT = -4.8


import pybamm
import numpy as np
import consts as c
from cell import Cell

model = pybamm.BaseModel()
geo = {}
cells = [Cell(f"Cell {i + 1}", model, geo) for i in range(NUM_CELLS)]

i_total = pybamm.Parameter("Input Current / Area") 

# Vcell1 - Vcell2 = 0
# (pos_phi1 - neg_phi1) - (pos_phi2 - neg_phi2) = 0
for i in range(len(cells) - 1):
    model.algebraic[cells[i].pos.phi] = cells[i + 1].pos.phi - cells[i].pos.phi
    model.algebraic[cells[i].neg.phi] = cells[i + 1].neg.phi - cells[i].neg.phi

    model.algebraic[cells[i].iapp] = cells[i].pos.phi_expr - cells[i].neg.phi_expr - (cells[i].pos.phi - cells[i].neg.phi)

    # temp3[cells[i].voltage] = cells[i + 1].voltage - cells[i].voltage
    # temp3[cells[i].iapp] = cells[i].pos.phi - cells[i].neg.phi - cells[i].voltage

model.algebraic[cells[-1].pos.phi] = cells[-1].pos.phi_expr - cells[-1].pos.phi
model.algebraic[cells[-1].neg.phi] = cells[-1].neg.phi_expr - cells[-1].neg.phi
model.algebraic[cells[-1].iapp] = i_total - sum(cell.iapp for cell in cells)


# best guesses
## current = input() / n
## pos_phi = p_OCP() @ t=0
## neg_phi = n_OCP() @ t=0
pos_phi_init = c.POS_OPEN_CIRCUIT_POTENTIAL(c.POS_CSN_INITIAL / c.POS_CSN_MAX)
neg_phi_init = c.NEG_OPEN_CIRCUIT_POTENTIAL(c.NEG_CSN_INITIAL / c.NEG_CSN_MAX)

parameters = {
    i_total.name: DISCHARGE_CURRENT
}

model.initial_conditions.update({
    **{ cell.pos.phi: pos_phi_init for cell in cells }, 
    **{ cell.neg.phi: neg_phi_init for cell in cells }, 
    **{ cell.iapp: DISCHARGE_CURRENT / NUM_CELLS for cell in cells }
})

for cell in cells:
    cell.set_parameters(parameters)

param_ob = pybamm.ParameterValues(parameters)
param_ob.process_model(model)
param_ob.process_geometry(geo)

PTS = 30
particles = [] # (pos1, pos2, pos3, pos4, neg1, neg2, neg3, neg4)
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

# solution.plot(list(model.variables.keys()))

### SOLVED
from matplotlib import pyplot as plt

fig, axs = plt.subplots(1, 4, figsize=(13,13))

for i in range(len(axs)):
    ax = axs[i]
    iapp = solution[cells[i].iapp.name].entries
    ax.plot(solution.t, np.round(iapp, 5))
    ax.set_xlabel("t")
    ax.set_ylabel("I (A)")
    ax.set_title(cells[i].iapp.name)

plt.tight_layout()
plt.show()