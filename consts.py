# assuming "Diffusion Coefficient" is constant w/ respect to Concentration @ r
D = 3.9e-14 # pybamm.Parameter("pDiffusion Coefficient")
F = 96485

# not scaling radius at the moment... 
R = 5.5e-06 

POS_CSN_MAX = 63104.0
POS_CSN_INITIAL = 2.5e4
POS_ELEC_THICKNESS = 7.56e-05
POS_ELEC_POROSITY = 0.33

NEG_CSN_MIN = 2000
NEG_CSN_INITIAL = 1.5e4
NEG_ELEC_THICKNESS = 8.52e-05
NEG_ELEC_POROSITY = 0.23

C_RATE = 5
RUNTIME_HOURS = 1 / C_RATE # hours

# RUNTIME_HOURS = 5