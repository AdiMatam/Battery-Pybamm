from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE
import pandas as pd
from matplotlib import pyplot as plt
from src.pack import Pack, Protocol, PROTOCOL_NAMES
import re


df = pd.read_csv("data/test_10/data.csv", index_col=[0,1])
cols = [col for col in list(df.columns) if "Anode" in col and "1,1" in col and "Current" in col]

## hardcoded

# side, jval
cols = cols[1:]

CYCLE = df.index.get_level_values(0)
PROTOCOL = df.index.get_level_values(1)

cc1 = df.loc[PROTOCOL.str.contains("CC_Charge")]
cv1 = df.loc[PROTOCOL.str.contains("CV_Charge")]

fig, ax = plt.subplots(2)

for c in [1]:
    cc = cc1.loc[c]
    cv = cv1.loc[c]

    for i in range(len(cols)):
        col = cols[i]
        ax[i].plot(cc['Clock'], cc[col], label=f"CC_{c}-{col}")
        ax[i].plot(cv['Clock'], cv[col], label=f"CV_{c}-{col}")
        ax[i].grid(True)        # turn grid on for this Axes
        ax[i].legend()          # draw this Axes’ legend

    #end = cv['Time'].tolist()[-1] + 10
    # plt.axvline(x=end, color='black', linestyle='--')

plt.show()
    