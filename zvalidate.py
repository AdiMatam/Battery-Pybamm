import pandas as pd
from matplotlib import pyplot as plt
import params as p

seidf = pd.read_csv("ANODE_2.csv")
#catdf = pd.read_csv("CATHODE_SEI_2.csv")

fig, ax1 = plt.subplots()

ax1.set_xlabel('Time (s)')
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# color = 'tab:red'
# ax1.set_ylabel('Current (A/m2)', color=color)

a = (3 * (1 - p.NEG_ELEC_POROSITY.get_value())) / p.PARTICLE_RADIUS.get_value()

#ax1.plot(seidf['Time'], seidf['I_INT'] * a * p.NEG_ELEC_THICKNESS.get_value(), color='blue', label='Intercalation')
#ax1.plot(seidf['Time'], seidf['I_INT'], color='blue', label='Intercalation')
#ax1.plot(seidf['Time'], seidf['I_SEI'], color='red', label='SEI')
#ax1.tick_params(axis='y', labelcolor='black')

# ax1.set_ylabel('Length (m)', color='black')
# ax1.plot(seidf['Time'], seidf['Length'], color='black')
# ax1.tick_params(axis='y', labelcolor='black')

#print(seidf['I_SEI'])
ax1.set_ylabel('Concentration', color='green')
ax1.plot(seidf['Time'], seidf['Anode vConcentration'], color='green')
ax1.tick_params(axis='y', labelcolor='green')

ax2.set_ylabel('Voltage (V)', color='tab:blue')  # we already handled the x-label with ax1
ax2.plot(seidf['Time'], seidf['Anode vPhi'], color='green')
# ax2.plot(seidf['Time'], catdf['Voltage'] - seidf['Voltage'], color='blue')
ax2.tick_params(axis='y', labelcolor='tab:blue')

plt.grid(linewidth=0.2)
plt.tight_layout()
plt.legend()
plt.show()