import random
from marquis import lico2_electrolyte_exchange_current_density_Dualfoil1998 as j0p
from marquis import graphite_electrolyte_exchange_current_density_Dualfoil1998 as j0n
from marquis import lico2_ocp_Dualfoil1998 as Up
from marquis import graphite_mcmb2528_ocp_Dualfoil1998 as Un

class VariatedParameter:
    def __init__(self, value: float, low_high: tuple):
        self.value = value
        self.low_high = low_high

    @classmethod
    def from_percent(cls, value: float, percent: float=1):
        offset = value * (percent / 100)
        return cls(value, (value - offset, value + offset))

    def rand_sample(self):
        return random.uniform(self.low_high[0], self.low_high[1]) 

    def get_value(self):
        return self.value











POS_CSN_MAX         = VariatedParameter.from_percent(51218, 0)  
POS_CSN_INITIAL     = VariatedParameter.from_percent(30730, 0)
POS_ELEC_THICKNESS  = VariatedParameter.from_percent(0.0001, 0)
POS_ELEC_POROSITY   = VariatedParameter.from_percent(0.50, 0)   

NEG_CSN_MIN         = VariatedParameter.from_percent(5027, 0)
NEG_CSN_INITIAL     = VariatedParameter.from_percent(19986, 0)
NEG_CSN_MAX         = VariatedParameter.from_percent(24983, 0)   
NEG_ELEC_THICKNESS  = VariatedParameter.from_percent(0.0001, 0)
NEG_ELEC_POROSITY   = VariatedParameter.from_percent(0.40, 0)

ELECTROLYTE_CONC    = VariatedParameter.from_percent(1000, 0)

PARTICLE_RADIUS     = VariatedParameter.from_percent(1e-05, 0)

POS_J0 = j0p
POS_OCP = Up

NEG_EXCHANGE_CURRENT_DENSITY = j0n
NEG_OCP = Un
