# assuming "Diffusion Coefficient" is constant w/ respect to Concentration @ r
# D = 2.0e-14 #3.9e-14
# D = 3.9e-14 # pybamm.Parameter("pDiffusion Coefficient")
F = 96485

# not scaling radius at the moment... 
R_GAS = 8.314
T = 298 # kelvin
RTF = R_GAS * T / F

def BIND_VALUES(fulldict: dict, subdict: dict):
    # absorption of particle parameters
    fulldict.update(
        {key.name : value for key, value in subdict.items()}
    )

def SET_MODEL_VARS(model, variables: list):
    model.variables.update(
        {var.name: var for var in variables}
    )

def SET_OUTPUTS(outputs: list, variables: list):
    outputs.extend( [var.name for var in variables] )