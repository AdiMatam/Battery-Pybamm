C_RATE = 1/20 
RUNTIME_HOURS = 1 / C_RATE # hours
NUM_CELLS = 1
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


#--------------------------------------------------

import consts as c


time_steps = len(solution.t)

voltages = []
for i in range(time_steps):
    voltages.append(0)

data = {
    'J0p': [],
    'up': [],
    'J0n': [],
    'un': []
}


for cell in cells:

    # length of entries == # of time steps (600)
    # surface concentration @ each time step
    surf_p = solution[cell.pos.surf_csn_name].entries[-1]
    surf_n = solution[cell.neg.surf_csn_name].entries[-1]

    ## j (electrode current density is constant throughout?)
    j_p = (1 * DISCHARGE_CURRENT) / (cell.pos_elec_thickness * cell.pos_elec_porosity) # solution[positive.j_name].entries[0]
    j_n = (-1 * DISCHARGE_CURRENT) / (cell.neg_elec_thickness * cell.neg_elec_porosity) # solution[negative.j_name].entries[0]

    for i in range(time_steps):
        # get surface concentration @ each time step
        inst_surf_p = surf_p[i]
        scaled_surf_p = inst_surf_p / cell.pos_csn_maxval
        # get current j0 (at i-th timestep)
        j0_p = p.POS_J0(p.ELECTROLYTE_CONC.get_value(), inst_surf_p, cell.pos_csn_maxval)

        # jp is negative. overpotential decreases when pos electrode being lithiated, 
        # so correct signs
        volmer_p = 2 * c.RTF * pybamm.arcsinh(j_p / (2 * j0_p)) # im getting higher overpot.
        
        # default function (given in pybamm basicSPM code -- check Up.py)
        up = p.Up(scaled_surf_p)

        data['J0p'].append(j0_p)
        data['up'].append(up)
        
        # -------------------------------------
        
        inst_surf_n = surf_n[i]
        scaled_surf_n = inst_surf_n / cell.neg_csn_maxval
        j0_n = p.NEG_J0(p.ELECTROLYTE_CONC.get_value(), inst_surf_n, cell.neg_csn_maxval)
        
        # jn is positive. overpotential increases when neg electrode being de-lithiated (discharge),
        # so, correct signs.
        volmer_n = 2 * c.RTF * pybamm.arcsinh(j_n / (2 * j0_n))

        # default function (given in pybamm basicSPM code -- check Un.py)
        un = p.Un(scaled_surf_n)
        print(un)

        data['J0n'].append(j0_n)
        data['un'].append(un)
        
        pos_v = up + volmer_p 
        v = pos_v + volmer_n - un
        voltages[i] += v.value # going from pybamm.Scalar() object to normal integer

import pandas as pd

df = pd.DataFrame(data)
print(df)
df.to_csv("../hopeful.csv", index=False)
