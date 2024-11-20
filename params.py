from math import exp
from src.variator import Variator

## .from_percent(value, % variation)
POS_DIFFUSION       = Variator.from_percent("Cathode Diffusion",       1.0e-14,     0)
POS_CSN_MAX         = Variator.from_percent("Cathode Max Concentration", 51555,     0)  
POS_CSN_INITIAL     = Variator.from_percent("Cathode Initial SOC",      51555*0.5,  0)
POS_ELEC_THICKNESS  = Variator.from_percent("Cathode Thickness",        80e-6,      0)
POS_ELEC_POROSITY   = Variator.from_gaussian_stddev("Cathode Porosity", 0.385,      0.0)   

NEG_DIFFUSION       = Variator.from_percent("Anode Diffusion",          2.0e-14,    0)
NEG_CSN_MAX         = Variator.from_percent("Anode Max Concentration",  30555,      0)   
NEG_CSN_INITIAL     = Variator.from_percent("Anode Initial SOC",       30555*0.74,  0)
NEG_ELEC_THICKNESS  = Variator.from_percent("Anode Thickness",          88e-6,      0)
NEG_ELEC_POROSITY   = Variator.from_gaussian_stddev("Anode Porosity",   0.485,      0.0)

PARTICLE_RADIUS     = Variator.from_percent("Particle Radius",          2e-06,      0)









import pybamm

def NEG_OCP(sto):
    x = sto
    return 0.7222 + 0.1387*x + 0.029*x**0.5 - 0.0172/x + 0.0019/(x**1.5) + 0.2808*exp(0.9-15*x) - 0.7984*exp(0.4465*x - 0.4108)

def NEG_OCP2(sto):
    x = sto
    return 0.7222 + 0.1387*x + 0.029*x**0.5 - 0.0172/x + 0.0019/(x**1.5) + 0.2808*pybamm.exp(0.9-15*x) - 0.7984*pybamm.exp(0.4465*x - 0.4108)

def POS_OCP(sto):
    y = sto
    num = -4.656 + 88.669*y**2 - 401.119*y**4 + 342.909*y**6 - 462.471*y**8 + 433.434*y**10
    den = -1 + 18.933*y**2 - 79.532*y**4 + 37.311*y**6 - 73.083*y**8 + 95.96*y**10
    return num / den


if __name__ == '__main__':
    for _ in range(10):
        print(POS_ELEC_POROSITY.sample())