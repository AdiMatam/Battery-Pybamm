NUM_CELLS = 2
NUM_CYCLES = 1
base_current = 1.20276592916666664
I_TOTAL = base_current * NUM_CELLS
VOLTAGE_CUTOFF = (2.0, 5.0) # effectively disabled

DISCRETE_PTS = 50
TIME_PTS = 250
RUNTIME_HOURS = 20






import pybamm
import params as p
from pack import Pack

i_param = pybamm.Parameter("Input Current / Area") 

model = pybamm.BaseModel()
geo = {}
parameters = {i_param.name: "[input]"}

pack = Pack(NUM_CELLS, model, geo, parameters, i_param, VOLTAGE_CUTOFF)

for i in range(2):
    for j in range(2):
        model.initial_conditions.update({
            pack.cells[i, j].iapp: -I_TOTAL / 2
        })

pack.build(DISCRETE_PTS)

variables = []
for cell in pack.flat_cells:
    variables.extend([
        cell.pos.surf_csn_name, 
        cell.neg.surf_csn_name, 
        # cell.voltage_name, 
        cell.iapp_name,
    ])

df = pack.cycler(I_TOTAL, NUM_CYCLES, RUNTIME_HOURS, TIME_PTS, variables, output_path="full_cycle_data.csv")

# from plotter import plot
# plot(df, cells)