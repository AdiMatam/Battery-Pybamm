# assuming "Diffusion Coefficient" is constant w/ respect to Concentration @ r
D = 3.9e-14 # pybamm.Parameter("pDiffusion Coefficient")
F = 96485

# not scaling radius at the moment... 
R = 5.22e-06 
R_GAS = 8.314
T = 298 # kelvin

POS_CSN_MAX = 40104.0
POS_CSN_INITIAL = 17038

POS_ELEC_THICKNESS = 7.56e-05
POS_ELEC_POROSITY = 0.33

#-------------------

NEG_CSN_MIN = 5000
NEG_CSN_INITIAL = 29866
NEG_CSN_MAX = 33133.0

NEG_ELEC_THICKNESS = 8.52e-05
NEG_ELEC_POROSITY = 0.25

C_RATE = 1
RUNTIME_HOURS = 1 / C_RATE # hours

# RUNTIME_HOURS = 10