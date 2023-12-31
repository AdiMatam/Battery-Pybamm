from marquis import lico2_electrolyte_exchange_current_density_Dualfoil1998 as j0p
from marquis import graphite_electrolyte_exchange_current_density_Dualfoil1998 as j0n
from marquis import lico2_ocp_Dualfoil1998 as Up
from marquis import graphite_mcmb2528_ocp_Dualfoil1998 as Un

# assuming "Diffusion Coefficient" is constant w/ respect to Concentration @ r
D = 3.9e-14 # pybamm.Parameter("pDiffusion Coefficient")
F = 96485

# not scaling radius at the moment... 
R = 1e-05
R_GAS = 8.314
T = 298 # kelvin

POS_CSN_MAX = 51218  # 40104.0
POS_CSN_INITIAL = 30730

POS_ELEC_THICKNESS = 0.0001
POS_ELEC_POROSITY = 0.50
POS_EXCHANGE_CURRENT_DENSITY = j0p
POS_OPEN_CIRCUIT_POTENTIAL = Up

#-------------------

NEG_CSN_MIN = 5026.70174550953 # 4354.87 #2500 # 4354.87 #2500
NEG_CSN_INITIAL = 19986
NEG_CSN_MAX = 24983

NEG_ELEC_THICKNESS = 0.0001
NEG_ELEC_POROSITY = 0.40
 
NEG_EXCHANGE_CURRENT_DENSITY = j0n
NEG_OPEN_CIRCUIT_POTENTIAL = Un

#-------------------
ELECTROLYTE_CONC = 1000

C_RATE = 1 
RUNTIME_HOURS = 1 / C_RATE # hours

# RUNTIME_HOURS = 10