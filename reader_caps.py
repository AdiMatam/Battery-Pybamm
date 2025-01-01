import pandas as pd
from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE
from matplotlib import pyplot as plt

## LOADING AN EXPERIMENT (folder name within data/)
noaged = Experiment("../data/S_1.0C_100_noage")
noaged.select_cycles(
    cycles=[1] + list(range(10,100,10))
)
noaged.select_attributes(["Cathode.*?Concentration"])

### Plot the CURRENT dataset (i.e. after all predecessing filters)
### isolate_cycles =True:  Plot data for EACH cycle as separate line with respect to "local time"
###                =False: Plot data with respect to "global time" (no delineation by cycle #)
noaged.plotter(isolate_cycles=True)

plt.legend()

plt.grid(True)
plt.show()