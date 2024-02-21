NUM_CELLS = 1
NUM_CYCLES = 1
base_current = 1.20276592916666664
I_TOTAL = base_current * NUM_CELLS
VOLTAGE_CUTOFF = (3.0, 5.0)

DISCRETE_PTS = 50
TIME_PTS = 250
RUNTIME_HOURS = 20






import pybamm
import params as p
from pack import Pack

i_total = pybamm.Parameter("Input Current / Area") 

model = pybamm.BaseModel()
geo = {}
parameters = {i_total.name: "[input]"}

pack = Pack(NUM_CELLS, model, geo, parameters, i_total, VOLTAGE_CUTOFF)
cells = pack.get_cells()

## pos_phi = p_OCP() @ t=0
## neg_phi = n_OCP() @ t=0
model.initial_conditions.update({
    **{ cell.pos.phi: p.POS_OCP(cell.pos_csn_ival / cell.pos_csn_maxval) for cell in cells }, 
    **{ cell.neg.phi: p.NEG_OCP(cell.neg_csn_ival / cell.neg_csn_maxval) for cell in cells }, 
    **{ cell.iapp:    -I_TOTAL / pack.get_num_cells() for cell in cells }
})

pack.build(DISCRETE_PTS)

variables = []
for cell in cells:
    variables.extend([
        cell.pos.surf_csn_name, 
        cell.neg.surf_csn_name, 
        cell.voltage_name, 
        cell.iapp_name
    ])

df = pack.cycler(I_TOTAL, NUM_CYCLES, RUNTIME_HOURS, TIME_PTS, variables, output_path="full_cycle_data.csv")

from plotter import plot
plot(df, cells)