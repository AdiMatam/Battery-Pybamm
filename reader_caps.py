import pandas as pd
from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE
from matplotlib import pyplot as plt

## LOADING AN EXPERIMENT (folder name within data/)
cc = Experiment("../data/100_noage_cc")
cc.select_cycles(
    cycles=[10, 30, 50]
)
cc.select_attributes(["Capacity"])

single = Experiment("../data/100_noage_cc_SINGLEMODE")
single.select_cycles(
    cycles=[10, 30, 50]
)
single.select_attributes(["Capacity"])

cc.plotter(isolate_cycles=True)
single.plotter(isolate_cycles=True)

plt.legend()
plt.grid(True)
plt.show()