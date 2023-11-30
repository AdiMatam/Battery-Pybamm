import pybamm
import numpy as np
import consts as c
from single_particle import SingleParticle

model = pybamm.BaseModel()

name = "Cell1"
pos1_iapp = pybamm.Variable(name + " Iapp")
pos1 = SingleParticle(name + " Pos Particle", +1, pos1_iapp)
neg1 = SingleParticle(name + " Neg Particle", -1, pos1_iapp)
pos_phi1 = pybamm.Variable(name + " Pos Phi")
neg_phi1 = pybamm.Variable(name + " Neg Phi")
pos1.process_model(model)
neg1.process_model(model)

name = "Cell2"
pos2_iapp = pybamm.Variable(name + " Iapp")
pos2 = SingleParticle(name + " Pos Particle", +1, pos2_iapp)
neg2 = SingleParticle(name + " Neg Particle", -1, pos2_iapp)
pos_phi2 = pybamm.Variable(name + " Pos Phi")
neg_phi2 = pybamm.Variable(name + " Neg Phi")
pos2.process_model(model)
neg2.process_model(model)

volt_force = pybamm.Variable("Enforced Voltage")
i_total = pybamm.Parameter("Input Current / Area") 

RTF = c.R_GAS * c.T / c.F
model.algebraic = {
    pos_phi1: pos1.u_func(pos1.surf_conc / pos1.conc_max) + (2*RTF*pybamm.arcsinh(pos1.j / (2 * pos1.j0))) - pos_phi1,
    neg_phi1: neg1.u_func(neg1.surf_conc / neg1.conc_max) + (2*RTF*pybamm.arcsinh(neg1.j / (2 * neg1.j0))) - neg_phi1,

    pos_phi2: pos2.u_func(pos2.surf_conc / pos2.conc_max) + (2*RTF*pybamm.arcsinh(pos2.j / (2 * pos2.j0))) - pos_phi2,
    neg_phi2: neg2.u_func(neg2.surf_conc / neg2.conc_max) + (2*RTF*pybamm.arcsinh(neg2.j / (2 * neg2.j0))) - neg_phi2,

    volt_force: (pos_phi1 - neg_phi1) - (pos_phi2 - neg_phi2),
    pos2_iapp: i_total - pos1_iapp
}