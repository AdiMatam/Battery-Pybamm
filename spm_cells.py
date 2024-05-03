NUM_CELLS = 2
NUM_CYCLES = 1
BASE_CURRENT = 1.20276592916666664
I_TOTAL = BASE_CURRENT * NUM_CELLS
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

string1 = pybamm.Variable("String 1 Iapp")
string2 = pybamm.Variable("String 2 Iapp")

c1volt = pybamm.Variable("Cell 1 Voltage")
c2volt = pybamm.Variable("Cell 2 Voltage")
c3volt = pybamm.Variable("Cell 3 Voltage")
c4volt = pybamm.Variable("Cell 4 Voltage")

c1 = Cell("C1", model, geo, parameters, string1, VOLTAGE_CUTOFF)
c2 = Cell("C2", model, geo, parameters, string1, VOLTAGE_CUTOFF)
c3 = Cell("C3", model, geo, parameters, string2, VOLTAGE_CUTOFF)
c4 = Cell("C4", model, geo, parameters, string2, VOLTAGE_CUTOFF)


ivp = lambda c: (p.POS_OCP(c.pos_csn_ival / c.pos_csn_maxval))
ivn = lambda c: (p.NEG_OCP(c.neg_csn_ival / c.neg_csn_maxval))
model.initial_conditions.update({
    string1: -BASE_CURRENT,
    string2: -BASE_CURRENT,
    c1volt: ivp(c1) - ivn(c1),
    c2volt: ivp(c2) - ivn(c2),
    c3volt: ivp(c3) - ivn(c3),
    c4volt: ivp(c4) - ivn(c4)
})

model.algebraic.update({
    string1: i_t - string2 - (string1),
    string2: c3volt + c4volt - (c1volt + c2volt),
    c1volt: c1.pos.bv - c1.neg.bv - c1volt, # 0 = pos_phi - neg_phi - c1volt
    c2volt: c2.pos.bv - c2.neg.bv - c2volt,
    c3volt: c3.pos.bv - c3.neg.bv - c3volt,
    c4volt: c4.pos.bv - c4.neg.bv - c4volt,
})

model.variables.update({
    string1.name: string1,
    string2.name: string2,
    c1volt.name: c1volt, # 0 = pos_phi - neg_phi - c1volt
    c2volt.name: c2volt,
    c3volt.name: c3volt,
    c4volt.name: c4volt,
})

for key in model.variables.keys():
    print(key)

particles = []
for c in (c1, c2, c3, c4):
    particles.append(c.pos)
    particles.append(c.neg)

param_ob = pybamm.ParameterValues(parameters)
param_ob.process_model(model)
param_ob.process_geometry(geo)

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
for t in time_steps:
    # print(solution[c1.pos.surf_csn_name](t)[-1])
    print(solution[string2.name](t))




# def iroot(c: Cell):
    # term =  2*R_GAS*T*pybamm.arcsinh(c.iapp / (2*c.neg.L*c.neg.a_term*c.neg.j0))
    # term -= 2*R_GAS*T*pybamm.arcsinh(c.iapp / (2*c.pos.L*c.pos.a_term*c.pos.j0))
    # term += (c.neg.ocp*F - c.pos.ocp*F + c.pote*F)

    # return term


# pack = Pack(NUM_CELLS, model, geo, parameters, i_param, VOLTAGE_CUTOFF)

# for i in range(2):
    # for j in range(2):
        # model.initial_conditions.update({
            # pack.cells[i, j].iapp: -I_TOTAL / 2
        # })

# pack.build(DISCRETE_PTS)

# variables = []
# for cell in pack.flat_cells:
    # variables.extend([
        # cell.pos.surf_csn_name, 
        # cell.neg.surf_csn_name, 
        # # cell.voltage_name, 
        # cell.iapp_name,
    # ])

# df = pack.cycler(I_TOTAL, NUM_CYCLES, RUNTIME_HOURS, TIME_PTS, variables, output_path="full_cycle_data.csv")

# from plotter import plot
# plot(df, cells)