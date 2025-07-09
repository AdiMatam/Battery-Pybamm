NUM_SERIES = 2
NUM_PARALLEL =2

USE_C_RATE = True
C_RATE = 1.0
TIME_PTS = 100
DISCRETE_PTS = 30
EXPERIMENT = "cccv_noage_200"


import pybamm
from src.pack import Pack, Protocol
pybamm.set_logging_level("WARNING")

model = pybamm.BaseModel()
geo = {}
parameters = {}

pack = Pack(EXPERIMENT, NUM_PARALLEL, NUM_SERIES, model, geo, parameters, aging=False)
pack.build(DISCRETE_PTS)

for i in range(200):
      pack.simulate(Protocol.CC_Discharge, 4000.0, c_rate=1.0, until=5.0)
      pack.simulate(Protocol.CC_Charge, 4000.0, c_rate=1.0, until=4.2*2)
      pack.simulate(Protocol.CV_Charge, 5000.0, until=27.2638366181154*2*0.1)
      pack.next_cycle()

pack.export_profile()