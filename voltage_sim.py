import pickle
import pybamm
import consts as c
from math import sqrt
from single_particle import SingleParticle

fle = open("Up_func.pkl", 'rb')
Up = pickle.load(fle)
fle.close()

fle = open("Un_func.pkl", 'rb')
Un = pickle.load(fle)
fle.close()

def post_process_voltage(solution: pybamm.Solution, positive: SingleParticle, negative: SingleParticle):
    global Up
    global Un

    ### NAIVE METHOD OF VOLTAGE CALCULATION (post processing)
    # Voltage = UP + VOL_P - UN - VOL_N

    voltages = []
    RTF = c.R_GAS * c.T / c.F

    surf_p = solution[positive.surf_conc_name].entries # length 24 (one at each time)
    surf_n = solution[negative.surf_conc_name].entries

    ## j (electrode current density is constant throughout?)
    j_p = solution[positive.j_name].entries[0]
    j_n = solution[negative.j_name].entries[0]

    time_steps = len(solution.t)

    for i in range(time_steps):
        inst_surf_p = surf_p[i]
        scaled_surf_p = inst_surf_p / c.POS_CSN_MAX
        j0_p = solution[positive.j0_name].entries[i] # sqrt(scaled_surf_p) * sqrt(1 - scaled_surf_p)

        volmer_p = 2 * RTF * pybamm.arcsinh(j_p / (2 * j0_p))
        up = Up(scaled_surf_p)
        
        inst_surf_n = surf_n[i]
        scaled_surf_n = inst_surf_n / c.NEG_CSN_MAX
        j0_n = solution[negative.j0_name].entries[i] # sqrt(scaled_surf_n) * sqrt(1 - scaled_surf_n)
        
        volmer_n = 2 * RTF * pybamm.arcsinh(j_n / (2 * j0_n))
        un = Un(scaled_surf_n)
        
        v = up - volmer_p - volmer_n - un # evaluation of pybamm.Scalar() objects
        voltages.append(v.value)

    return voltages
