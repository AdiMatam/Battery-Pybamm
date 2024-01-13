import pybamm
import numpy as np
import consts as c
from cell import Cell

model = pybamm.BaseModel()
geo = {}

cells = []
for i in range(1, 5):
    cells.append(Cell(f"Cell{i}", model, geo))

i_total = pybamm.Parameter("Input Current / Area") 

# Vcell1 - Vcell2 = 0
# (pos_phi1 - neg_phi1) - (pos_phi2 - neg_phi2) = 0
temp3 = {}
for i in range(len(cells) - 1):
    temp3[cells[i].voltage] = cells[i + 1].voltage - cells[i].voltage
    temp3[cells[i].iapp] = cells[i].pos.phi - cells[i].neg.phi - cells[i].voltage

temp3[cells[-1].voltage] = cells[-1].pos.phi - cells[-1].neg.phi - cells[-1].voltage
temp3[cells[-1].iapp] = i_total - sum(cell.iapp for cell in cells)

model.algebraic = temp3

# best guesses
## current = input() / n
## pos_phi = p_OCP() @ t=0
## neg_phi = n_OCP() @ t=0
pos_phi_init = c.POS_OPEN_CIRCUIT_POTENTIAL(c.POS_CSN_INITIAL / c.POS_CSN_MAX)
neg_phi_init = c.NEG_OPEN_CIRCUIT_POTENTIAL(c.NEG_CSN_INITIAL / c.NEG_CSN_MAX)

model.initial_conditions.update({**{cell.voltage: pos_phi_init - neg_phi_init for cell in cells}, **{cell.iapp: -1.2 for cell in cells}})

param_dict = {
    i_total.name: -4.8
}

for cell in cells:
    cell.set_parameters(param_dict)

param_ob = pybamm.ParameterValues(param_dict)
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

solver = pybamm.CasadiSolver(mode="safe")
# solver = pybamm.ScipySolver()
time_steps = np.linspace(0, 3600 * 20, 250)
solution = solver.solve(model, time_steps)

solution.plot(list(model.variables.keys()))
