import pickle
import pybamm
import consts as cc
from particle_anode import Anode
from marquis import lico2_ocp_Dualfoil1998 as Up
from marquis import graphite_mcmb2528_ocp_Dualfoil1998 as Un

def post_process_voltage(solution: pybamm.Solution, positive: Anode, negative: Anode):
    voltages = []
    RTF = cc.R_GAS * cc.T / cc.F

    # length of entries == # of time steps (600)
    # surface concentration @ each time step
    surf_p = solution[positive.surf_csn_name].entries
    surf_n = solution[negative.surf_csn_name].entries

    ## j (electrode current density is constant throughout?)
    j_p = solution[positive.j_name].entries[0]
    j_n = solution[negative.j_name].entries[0]

    time_steps = len(solution.t)

    for i in range(time_steps):
        # get surface concentration @ each time step
        inst_surf_p = surf_p[i]
        scaled_surf_p = inst_surf_p / cc.POS_CSN_MAX
        # get current j0 (at i-th timestep)
        j0_p = solution[positive.j0_name].entries[i]

        # jp is negative. overpotential decreases when pos electrode being lithiated, 
        # so correct signs
        volmer_p = 2 * RTF * pybamm.arcsinh(j_p / (2 * j0_p))
        
        # default function (given in pybamm basicSPM code -- check Up.py)
        up = Up(scaled_surf_p)
        
        # -------------------------------------
        
        inst_surf_n = surf_n[i]
        scaled_surf_n = inst_surf_n / cc.NEG_CSN_MAX
        j0_n = solution[negative.j0_name].entries[i]
        
        # jn is positive. overpotential increases when neg electrode being de-lithiated (discharge),
        # so, correct signs.
        volmer_n = 2 * RTF * pybamm.arcsinh(j_n / (2 * j0_n))

        # default function (given in pybamm basicSPM code -- check Un.py)
        un = Un(scaled_surf_n)
        
        pos_v = up + volmer_p 
        v = pos_v - (volmer_n + un)
        voltages.append(v.value) # going from pybamm.Scalar() object to normal integer

    return voltages
