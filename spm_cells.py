C_RATE = 1/20 
I_TOTAL = 0
NUM_CELLS = 4
DISCRETE_PTS = 30
TIME_PTS = 250

RUNTIME_HOURS = 1 / C_RATE # hours


import pybamm
import numpy as np
from cell import Cell

model = pybamm.BaseModel()
geo = {}

i_total = pybamm.Parameter("Input Current / Area") 

parameters = {}
cells = [Cell(f"Cell {i + 1}", model, geo, parameters) for i in range(NUM_CELLS)]

# Vcell1 - Vcell2 = 0
# (pos_phi1 - neg_phi1) - (pos_phi2 - neg_phi2) = 0
for i in range(len(cells) - 1):
    model.algebraic[cells[i].pos.phi] = cells[i + 1].pos.phi - cells[i].pos.phi
    model.algebraic[cells[i].neg.phi] = cells[i + 1].neg.phi - cells[i].neg.phi

    model.algebraic[cells[i].iapp] = cells[i].pos.bv_term - cells[i].neg.bv_term - (cells[i].pos.phi - cells[i].neg.phi)

model.algebraic[cells[-1].pos.phi] = cells[-1].pos.bv_term - cells[-1].pos.phi
model.algebraic[cells[-1].neg.phi] = cells[-1].neg.bv_term - cells[-1].neg.phi
model.algebraic[cells[-1].iapp] = i_total - sum(cell.iapp for cell in cells)

## CAPACITY CALCULATION
capacity = NUM_CELLS * min(cell.capacity for cell in cells) # this is in Ah/m^2
i_total_val = capacity / RUNTIME_HOURS
print("Current:", -i_total_val)

parameters[i_total.name] = I_TOTAL or -i_total_val # -10 # -i_total_val


# best guesses
## current = input() / n
## pos_phi = p_OCP() @ t=0
## neg_phi = n_OCP() @ t=0
model.initial_conditions.update({
    **{ cell.pos.phi: cell.pos_phi_init for cell in cells }, 
    **{ cell.neg.phi: cell.neg_phi_init for cell in cells }, 
    **{ cell.iapp: -i_total_val / NUM_CELLS for cell in cells }
})


for cell in cells:
    print(cell.pos_elec_porosity)


param_ob = pybamm.ParameterValues(parameters)
param_ob.process_model(model)
param_ob.process_geometry(geo)

particles = [] 
for cell in cells:
    particles.append(cell.pos)
    particles.append(cell.neg)

mesh = pybamm.Mesh(geo, 
    { p.domain: pybamm.Uniform1DSubMesh for p in particles },
    { p.r: DISCRETE_PTS for p in particles }
)

disc = pybamm.Discretisation(mesh, 
    { p.domain: pybamm.FiniteVolume() for p in particles }
)
disc.process_model(model)

solver = pybamm.CasadiSolver()
time_steps = np.linspace(0, 3600 * RUNTIME_HOURS, TIME_PTS)
solution = solver.solve(model, time_steps)

solution.plot(list(model.variables.keys()))

# SOLVED
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