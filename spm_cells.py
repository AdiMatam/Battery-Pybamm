C_RATE = 1/20 
I_TOTAL = 0
NUM_CELLS = 4
DISCRETE_PTS = 30
TIME_PTS = 250

RUNTIME_HOURS = 1 / C_RATE # hours


import pybamm
import numpy as np
from cell import Cell
import params as p

model = pybamm.BaseModel()
geo = {}

i_total = pybamm.Parameter("Input Current / Area") 

parameters = {}
parameters[i_total.name] = "[input]" 
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
i_from_capacity = capacity / RUNTIME_HOURS
print("Current:", -i_from_capacity)


## pos_phi = p_OCP() @ t=0
## neg_phi = n_OCP() @ t=0
model.initial_conditions.update({
    **{ cell.pos.phi: p.POS_OCP(cell.pos_csn_ival / cell.pos_csn_maxval) for cell in cells }, 
    **{ cell.neg.phi: p.NEG_OCP(cell.neg_csn_ival / cell.neg_csn_maxval) for cell in cells }, 
    **{ cell.iapp: -i_from_capacity / NUM_CELLS for cell in cells }
})

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

cycles = 4
solver = pybamm.CasadiSolver()
time_steps = np.linspace(0, 3600 * RUNTIME_HOURS, TIME_PTS)
total_time_steps = np.linspace(0, 3600 * RUNTIME_HOURS * cycles, TIME_PTS * cycles)

inps = {
    i_total.name: -i_from_capacity,
    ## intial conditions for particle concentrations are set here!!
    ** {cell.pos.c_0.name: cell.pos_csn_ival for cell in cells},
    ** {cell.neg.c_0.name: cell.neg_csn_ival for cell in cells},
}


import pandas as pd

keys = []
for cell in cells:
    keys.extend([
        cell.pos.surf_csn_name, 
        cell.neg.surf_csn_name, 
        cell.voltage_name, 
        cell.iapp_name
    ])

subdfs = []

solution = None
for i in range(cycles):
    solution = solver.solve(model, time_steps, inputs=inps)

    subdf = pd.DataFrame(columns=keys)

    for key in keys:
        data = solution[key].entries
        if "Concentration" in key:
            data = data[-1]

        subdf[key] = data

    subdfs.append(subdf)

    pos_next = solution[cells[0].pos.surf_csn_name].entries[-1][-1]
    neg_next = solution[cells[0].neg.surf_csn_name].entries[-1][-1]

    inps[i_total.name] *= -1

    inps.update({
        ** {cell.pos.c_0.name: pos_next for cell in cells},
        ** {cell.neg.c_0.name: neg_next for cell in cells},
    })

df = pd.concat(subdfs, ignore_index=True)
df.insert(0, 'Time', total_time_steps)
df.to_csv("full_cycle_data.csv", index=False)

from  matplotlib import pyplot as plt

fig, ax1 = plt.subplots()

ax1.set_xlabel('Time (s)')

color = 'tab:red'
ax1.set_ylabel('Concentration (mol / m2)', color=color)
ax1.plot(df['Time'], df[cells[0].pos.surf_csn_name], color=color)
ax1.plot(df['Time'], df[cells[0].neg.surf_csn_name], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Potential (V)', color=color)  # we already handled the x-label with ax1
ax2.plot(df['Time'], df[cells[0].voltage_name], color=color)
ax2.tick_params(axis='y', labelcolor=color)

plt.tight_layout()
plt.show()
