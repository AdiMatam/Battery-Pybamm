import pandas as pd
from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE
from matplotlib import pyplot as plt

my_cycles = list(range(1,31))

def extract_caps(experiment):
    # experiment = Experiment("100_age_experimentcv_CONTEXTSW")
    experiment.select_cycles(
        cycles=my_cycles,
        protocols=[DISCHARGE]
    )
    experiment.select_attributes(["Capacity"])
    CELLS = experiment.pack.cells
    capdf = experiment.data[CELLS[0,0].capacity.name]

    caps = []
    for i in my_cycles:
        caps.append( capdf[i].iloc[-1] )

    return caps

e1 = Experiment("30_age_cccv_CONTEXTSW")
e2 = Experiment("30_noage_cccv_CONTEXTSW")
pycap1 = extract_caps(e1)
pycap2 = extract_caps(e2)

aged = [pycap1]
noaged = [pycap2]


f = open("caps_age.txt", 'r')
for line in f:
    if not line.strip():
        continue
    if (line.startswith("-")):
        aged.append([])
    else:
        aged[-1].append(float(line.strip()))   

f = open("caps_noage.txt", 'r')
for line in f:
    if not line.strip():
        continue
    if (line.startswith("-")):
        noaged.append([])
    else:
        noaged[-1].append(float(line.strip()))


plt.scatter(my_cycles, aged[0], label="Pybamm_SEI")
plt.scatter(my_cycles, noaged[0], label="Pybamm_NoAge")

plt.scatter(my_cycles, aged[1], label="Sundials_SEI")
plt.scatter(my_cycles, noaged[1], label="Sundials_NoAge")

plt.scatter(my_cycles, aged[2], label="Matlab_SEI")
plt.scatter(my_cycles, noaged[2], label="Matlab_NoAge")

plt.legend()
plt.show()
