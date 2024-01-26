# assuming "Diffusion Coefficient" is constant w/ respect to Concentration @ r
D = 3.9e-14 # pybamm.Parameter("pDiffusion Coefficient")
F = 96485

# not scaling radius at the moment... 
R = 1e-05
R_GAS = 8.314
T = 298 # kelvin
RTF = R_GAS * T / F

