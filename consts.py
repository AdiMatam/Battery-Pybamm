# assuming "Diffusion Coefficient" is constant w/ respect to Concentration @ r
# D = 2.0e-14 #3.9e-14
# D = 3.9e-14 # pybamm.Parameter("pDiffusion Coefficient")
from constant_parameter import ConstantParameter


F = 96485
R_GAS = 8.314
T = 298 # kelvin
RTF = R_GAS * T / F

def BIND_VALUES(fulldict: dict, subdict: dict):
    for key, value in subdict.items():
        if isinstance(key, ConstantParameter) and type(value) is not str:
            key.value = value ## for querying purposes, the parameter will also store its assigned value

        fulldict[key.name] = value ## format for ParameterValues object (processed by pybamm internal)


def SET_MODEL_VARS(model, variables: list):
    model.variables.update(
        {var.name: var for var in variables}
    )

def SET_OUTPUTS(outputs: list, variables: list):
    outputs.extend( [var.name for var in variables] )