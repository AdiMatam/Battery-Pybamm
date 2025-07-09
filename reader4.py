import pandas as pd
from matplotlib import pyplot as plt
from src.pack import Pack, Protocol, PROTOCOL_NAMES
import re


#fig, ax = plt.subplots(1)

i = 0
for name in ("cccv_noage_200","cccv_noage_200_prime"):

    df = pd.read_csv(f"data/{name}/data.csv", index_col=[0,1])
    cols = [col for col in list(df.columns) if "Capacity" in col]

    CYCLE = df.index.get_level_values(0)
    PROTOCOL = df.index.get_level_values(1)

    dch = df.loc[PROTOCOL.str.contains("Discharge")]

    cycles = [2] + list(range(10, 210, 10))
    col = cols[0]
    caps = []
    for c in cycles:
        caps.append(dch.loc[c][col].iloc[-1])

    plt.scatter(cycles, caps, label=f"{name}")

    plt.legend()
    #i += 1
    
plt.savefig(f"CCCV_identical_modes.png")
plt.show()