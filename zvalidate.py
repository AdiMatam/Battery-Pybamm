import pandas as pd
from matplotlib import pyplot as plt
import params as p

seidf = pd.read_csv("ANODE_SEI_3.csv")
noseidf =pd.read_csv("ANODE_NOSEI_1.csv") 
catdf = pd.read_csv("CATHODE_1.csv")


fig, ax1 = plt.subplots()

ax1.set_xlabel('Time (s)')
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:red'
ax1.set_ylabel('Current (A/m2)', color=color)

a = (3 * (1 - p.NEG_ELEC_POROSITY.get_value())) / p.PARTICLE_RADIUS.get_value()

#ax1.plot(seidf['Time'], seidf['I_INT'] * a * p.NEG_ELEC_THICKNESS.get_value(), color='blue', label='Intercalation')
ax1.plot(seidf['Time'], seidf['I_SEI'], color='red', label='SEI')
ax1.tick_params(axis='y', labelcolor=color)

# ax1.set_ylabel('Length (m)', color=color)
# ax1.plot(seidf['Time'], seidf['Length'], color='blue', label='Intercalation')
# ax1.tick_params(axis='y', labelcolor=color)

# ax1.set_ylabel('Concentration', color=color)
# ax1.plot(seidf['Time'], noseidf['Concentration'], color='blue', label='No SEI')
# ax1.plot(seidf['Time'], seidf['Concentration'], color='red', label='SEI')
# ax1.tick_params(axis='y', labelcolor=color)

# ax2.set_ylabel('Voltage (V)', color='tab:blue')  # we already handled the x-label with ax1
# ax2.plot(seidf['Time'], noseidf['Voltage'], color='blue', label='NO SEI')
# ax2.plot(seidf['Time'], seidf['Voltage'], color='red', label='SEI')
# ax2.tick_params(axis='y', labelcolor='tab:blue')

# ax2.set_ylabel('Voltage (V)', color='tab:blue')  # we already handled the x-label with ax1
# ax2.plot(catdf['Time'], catdf['Voltage'], color='blue', label='NO SEI')
# ax2.tick_params(axis='y', labelcolor='tab:blue')

# ax2.set_ylabel('Voltage (V)', color='tab:blue')  # we already handled the x-label with ax1
# ax2.plot(catdf['Time'], seidf['Voltage'], color='blue', label='SEI')
# ax2.plot(catdf['Time'], noseidf['Voltage'], color='red', label='NO_SEI')
# ax2.tick_params(axis='y', labelcolor='tab:blue')


plt.grid(linewidth=0.2)
plt.tight_layout()
plt.legend()
plt.show()