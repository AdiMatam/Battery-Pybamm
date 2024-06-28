from variatedparam import VariatedParameter

VariatedParameter.OVERRIDE_VARIATON = True

from marquis import lico2_electrolyte_exchange_current_density_Dualfoil1998 as j0p
from marquis import graphite_electrolyte_exchange_current_density_Dualfoil1998 as j0n
from marquis import lico2_ocp_Dualfoil1998 as Up
from marquis import graphite_mcmb2528_ocp_Dualfoil1998 as Un

## .from_percent(value, % variation)
POS_CSN_MAX         = VariatedParameter.from_percent(51218, 0)  
POS_CSN_INITIAL     = VariatedParameter.from_percent(30730, 0)
POS_ELEC_THICKNESS  = VariatedParameter.from_percent(0.0001, 0)
POS_ELEC_POROSITY   = VariatedParameter.from_percent(0.50, 5)   

NEG_CSN_MIN         = VariatedParameter.from_percent(2027, 0)
NEG_CSN_INITIAL     = VariatedParameter.from_percent(19986, 0)
NEG_CSN_MAX         = VariatedParameter.from_percent(24983, 0)   
NEG_ELEC_THICKNESS  = VariatedParameter.from_percent(0.0001, 0.0)
NEG_ELEC_POROSITY   = VariatedParameter.from_percent(0.40, 5)

ELECTROLYTE_CONC    = VariatedParameter.from_percent(1000, 0)

PARTICLE_RADIUS     = VariatedParameter.from_percent(1e-05, 10)

POS_J0 = j0p
POS_OCP = Up

NEG_J0 = j0n
NEG_OCP = Un
