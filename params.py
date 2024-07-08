from variatedparam import VariatedParameter

VariatedParameter.OVERRIDE_VARIATON = True

from ocp import neg_ocp, pos_ocp

## .from_percent(value, % variation)
POS_DIFFUSION       = VariatedParameter.from_percent(1.0e-14, 0)
POS_CSN_MAX         = VariatedParameter.from_percent(51555, 0)  
POS_CSN_INITIAL     = VariatedParameter.from_percent(51555*0.5, 0)
POS_ELEC_THICKNESS  = VariatedParameter.from_percent(80e-6, 0)
POS_ELEC_POROSITY   = VariatedParameter.from_percent(0.385, 5)   

NEG_DIFFUSION       = VariatedParameter.from_percent(2.0e-14, 0)
NEG_CSN_MAX         = VariatedParameter.from_percent(30555, 0)   
NEG_CSN_INITIAL     = VariatedParameter.from_percent(30555*0.74, 0)
NEG_ELEC_THICKNESS  = VariatedParameter.from_percent(88e-6, 0.0)
NEG_ELEC_POROSITY   = VariatedParameter.from_percent(0.485, 5)

#ELECTROLYTE_CONC    = VariatedParameter.from_percent(1000, 0)

PARTICLE_RADIUS     = VariatedParameter.from_percent(2e-06, 10)

POS_OCP = pos_ocp
NEG_OCP = neg_ocp
