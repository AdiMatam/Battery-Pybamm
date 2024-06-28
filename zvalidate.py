import pandas as pd
from matplotlib import pyplot as plt

#seidf = pd.read_csv("cycles_SEI_FULL_1.csv")
seidf = pd.read_csv("cycles_ANO_1.csv")
noseidf = pd.read_csv("cycles_NOSEI_1.csv")
catdf = pd.read_csv("cycles_CATH_1.csv")


fig, ax1 = plt.subplots()

ax1.set_xlabel('Time (s)')
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:red'
# ax1.set_ylabel('Current (I/m2)', color=color)
# ax1.plot(seidf['Time'], seidf['I_INT'], color='red', label='Intercalation')
# ax1.plot(seidf['Time'], seidf['I_SEI'], color='blue', label='SEI')
# ax1.tick_params(axis='y', labelcolor=color)

ax2.set_ylabel('Voltage (V)', color='tab:blue')  # we already handled the x-label with ax1
ax2.plot(seidf['Time'], catdf['Voltage'] - seidf['Voltage'], color='blue', label='SEI')
#ax2.plot(seidf['Time'], catdf['Voltage'], color='red', label='NO SEI')
ax2.tick_params(axis='y', labelcolor='tab:blue')

"""
catdf['Voltage'] - 
catdf['Voltage'] - 
"""

plt.grid(linewidth=0.2)
plt.tight_layout()
plt.legend()
plt.show()