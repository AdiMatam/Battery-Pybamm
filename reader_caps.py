import pandas as pd
from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE
from matplotlib import pyplot as plt

## LOADING AN EXPERIMENT (folder name within data/)
cc = Experiment("../data/100_noage_cccv_CONTEXTSW_CRAZY")
cc.select_cycles(
    cycles=[10, 50, 90]
)
cc.select_attributes(["Pack.*?Voltage"])

cc.plotter(isolate_cycles=True)

plt.legend()
plt.grid(True)
plt.show()