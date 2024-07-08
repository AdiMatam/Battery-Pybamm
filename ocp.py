import pybamm
from matplotlib import pyplot as plt
import numpy as np
import math

## wrong initial ocp?

def neg_ocp(sto):
    x = sto
    return 0.7222 + 0.1387*x + 0.029*x**0.5 - 0.0172/x + 0.0019/(x**1.5) + 0.2808*pybamm.exp(0.9-15*x) - 0.7984*pybamm.exp(0.4465*x - 0.4108)


def neg1_ocp(sto):
    x = sto
    return 0.7222 + 0.1387*x + 0.029*x**0.5 - 0.0172/x + 0.0019/(x**1.5) + 0.2808*math.exp(0.9-15*x) - 0.7984*math.exp(0.4465*x - 0.4108)

def pos_ocp(sto):
    y = sto
    num = -4.656 + 88.669*y**2 - 401.119*y**4 + 342.909*y**6 - 462.471*y**8 + 433.434*y**10
    den = -1 + 18.933*y**2 - 79.532*y**4 + 37.311*y**6 - 73.083*y**8 + 95.96*y**10
    return num / den

import params as p

if __name__ == '__main__':
    x = np.linspace(0.5, 1.0, 1000)
    y = np.linspace(0.75, 0.01, 1000)
    p_ocps = [pos_ocp(i) for i in x]
    n_ocps = [neg1_ocp(i) for i in y]

    #plt.plot(x, p_ocps)
    plt.plot(y, n_ocps)
    plt.show()

    # print(neg_ocp(p.NEG_CSN_INITIAL.get_value() / p.NEG_CSN_MAX.get_value()))
    # print(pos_ocp(p.POS_CSN_INITIAL.get_value() / p.POS_CSN_MAX.get_value()))
    
