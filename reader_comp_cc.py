import pandas as pd
from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE
from matplotlib import pyplot as plt

## LOADING AN EXPERIMENT (folder name within data/)
cc = Experiment("100_noage_cc_CONTEXTSW")
cc.select_cycles(
    cycles=[1],
    protocols=[DISCHARGE, CC_CHARGE]
)
cc.select_attributes(["Concentration"])
cc.plotter(isolate_cycles=True)

t = []
anode_conc = []
cathode_conc = []
voltage = []

with open("suncycle.txt") as f:
    for line in f:
        data = line.split("|")[:-1]
        t.append(float(data[0].strip()))
        cathode_conc.append(float(data[1].strip()))
        anode_conc.append(float(data[2].strip()))
        voltage.append(float(data[5].strip()))

plt.plot(t, cathode_conc, label="Sundials_Cathode")
plt.plot(t, anode_conc, label="Sundials_Anode")
#plt.plot(t, voltage, label="Sundials_Voltage")

plt.legend()
plt.grid(True)
plt.show()