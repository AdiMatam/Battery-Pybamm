import pandas as pd
from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE
from matplotlib import pyplot as plt

MY_CYCLES=[10, 50, 90, 100]

## LOADING AN EXPERIMENT (folder name within data/)
cc = Experiment("100_noage_cccv_CONTEXTSW")
cc.select_cycles(
    cycles=MY_CYCLES,
)
cc.select_attributes(["Capacity"])


t = []
anode_conc = []
cathode_conc = []
voltage = []
capacity = []

old_cycle = 0

with open("spm_dns.txt") as f:
    for line in f:
        if line.startswith("Cycle"):
            cycle = int(line.strip().split(" ")[-1])
            if old_cycle != cycle:
                t.append([])
                anode_conc.append([])
                cathode_conc.append([])
                voltage.append([])
                capacity.append([])
                old_cycle = cycle
            
        else:
            data = line.split("|")[:-1]
            t[old_cycle-1].append(float(data[0].strip()))
            cathode_conc[old_cycle-1].append(float(data[1].strip()))
            anode_conc[old_cycle-1].append(float(data[2].strip()))
            voltage[old_cycle-1].append(float(data[3].strip()))
            capacity[old_cycle-1].append(float(data[4].strip()))
            


fig, ax = plt.subplots(2)

# plt.plot(t, cathode_conc, label="Sundials_Cathode")
# plt.plot(t, anode_conc, label="Sundials_Anode")
for c in MY_CYCLES:
    plot_t = [t[c-1][i] - t[c-1][0] for i in range(len(t[c-1]))]
    ax[1].plot(plot_t, capacity[c-1], label=f"C{c}_Sundials")

# cc.plotter(isolate_cycles=True, plot=ax[0])

ax[0].legend()
ax[1].legend()
ax[0].grid(True)
ax[1].grid(True)
plt.show()