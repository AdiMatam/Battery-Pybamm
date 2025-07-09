from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE
import pandas as pd
from matplotlib import pyplot as plt
from src.pack import Pack, Protocol, PROTOCOL_NAMES
import re


df = pd.read_csv("data/test_10/data.csv", index_col=[0,1])
cols = [col for col in list(df.columns) if "SEI" in col and "1,1" in col]
col = cols[0]

CYCLE = df.index.get_level_values(0)
PROTOCOL = df.index.get_level_values(1)

#cyc = CYCLE.isin([1])
#cc = df.loc[cyc & PROTOCOL.str.contains("CC_Charge")]
#cv = df.loc[cyc & PROTOCOL.str.contains("CV_Charge")]

# for col in cols:
    # label = re.search(r"(Cell (\d+),(\d+)).*?", col).group()
    # plt.plot(cc['Clock'], cc[col], label=f"CC_{label}")
    # plt.plot(cv['Clock'], cv[col], label=f"CV_{label}")
    # print(cv['Clock'].tolist()[-1])

cc1 = df.loc[PROTOCOL.str.contains("CC_Charge")]
cv1 = df.loc[PROTOCOL.str.contains("CV_Charge")]

for c in [1, 2, 3]:
    cc = cc1.loc[c]
    cv = cv1.loc[c]

    plt.plot(cc['Clock'], cc[col], label=f"CC_{c}")
    plt.plot(cv['Clock'], cv[col], label=f"CV_{c}")

    end = cv['Clock'].tolist()[-1] + 10
    plt.axvline(x=end, color='black', linestyle='--')

plt.grid(True)
plt.legend()
plt.show()
    