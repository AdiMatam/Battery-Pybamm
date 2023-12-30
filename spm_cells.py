import pybamm
import numpy as np
import consts as c
from single_particle import SingleParticle

model = pybamm.BaseModel()

name = "Cell1"
cell1_iapp = pybamm.Variable(name + " Iapp")
pos1 = SingleParticle(name + " Pos Particle", +1, cell1_iapp)
neg1 = SingleParticle(name + " Neg Particle", -1, cell1_iapp)
pos_phi1 = pybamm.Variable(name + " Pos Phi")
neg_phi1 = pybamm.Variable(name + " Neg Phi")
pos1.process_model(model)
neg1.process_model(model)

name = "Cell2"
cell2_iapp = pybamm.Variable(name + " Iapp")
pos2 = SingleParticle(name + " Pos Particle", +1, cell2_iapp)
neg2 = SingleParticle(name + " Neg Particle", -1, cell2_iapp)
pos_phi2 = pybamm.Variable(name + " Pos Phi")
neg_phi2 = pybamm.Variable(name + " Neg Phi")
pos2.process_model(model)
neg2.process_model(model)

i_total = pybamm.Parameter("Input Current / Area") 

RTF = c.R_GAS * c.T / c.F

# Vcell1 - Vcell2 = 0
# (pos_phi1 - neg_phi1) - (pos_phi2 - neg_phi2) = 0
# write as: 
#       pos_phi1: (neg_phi1 + pos_phi2 - neg_phi2) - pos_phi1

model.algebraic = {
    pos_phi1: (neg_phi1 + pos_phi2 - neg_phi2) - pos_phi1,

    cell1_iapp: pos1.ocp + 2*RTF*pybamm.arcsinh((cell1_iapp / (pos1.L * pos1.a_term)) / (2 * pos1.j0)) - pos_phi1,

    neg_phi1: neg1.ocp + (2*RTF*pybamm.arcsinh(neg1.j / (2 * neg1.j0))) - neg_phi1,
    pos_phi2: pos2.ocp + (2*RTF*pybamm.arcsinh(pos2.j / (2 * pos2.j0))) - pos_phi2,
    neg_phi2: neg2.ocp + (2*RTF*pybamm.arcsinh(neg2.j / (2 * neg2.j0))) - neg_phi2,

    cell2_iapp: (i_total - cell1_iapp) - cell2_iapp
}

# best guesses
## current = input() / 2
## pos_phi = p_OCP() @ t=0
## neg_phi = n_OCP() @ t=0

model.initial_conditions.update({
    pos_phi1: 4.027031539154782,
    neg_phi1: 0.17519402534096734,

    pos_phi2: 4.027031539154782,
    neg_phi2: 0.17519402534096734,

    cell1_iapp: -1.2,
    cell2_iapp: -1.2,
})

model.variables.update({
    pos_phi1.name: pos_phi1,
    neg_phi1.name: neg_phi1,
    pos_phi2.name: pos_phi2,
    neg_phi2.name: neg_phi2,           

    cell1_iapp.name: cell1_iapp,
    cell2_iapp.name: cell2_iapp,
})

geo = {}
pos1.process_geometry(geo)
neg1.process_geometry(geo)

pos2.process_geometry(geo)
neg2.process_geometry(geo)

param_dict = {
    i_total.name: -2.4
}

for p in (pos1, pos2):
    p.process_parameters(param_dict, {
        p.conc_0:    c.POS_CSN_INITIAL,
        p.L:         c.POS_ELEC_THICKNESS,
        p.eps_n:     c.POS_ELEC_POROSITY,
        p.conc_max:  c.POS_CSN_MAX,
        p.j0:        c.POS_EXCHANGE_CURRENT_DENSITY,
        p.ocp:       c.POS_OPEN_CIRCUIT_POTENTIAL
    })

for n in (neg1, neg2):
    n.process_parameters(param_dict, {
        n.conc_0:    c.NEG_CSN_INITIAL,
        n.L:         c.NEG_ELEC_THICKNESS,
        n.eps_n:     c.NEG_ELEC_POROSITY,
        n.conc_max:  c.NEG_CSN_MAX,
        n.j0:        c.NEG_EXCHANGE_CURRENT_DENSITY,
        n.ocp:       c.NEG_OPEN_CIRCUIT_POTENTIAL
    })

param = pybamm.ParameterValues(param_dict)
param.process_model(model)
param.process_geometry(geo)

PTS = 30
particles = (pos1, pos2, neg1, neg2)
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