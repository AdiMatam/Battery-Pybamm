from tkinter import W
from  matplotlib import pyplot as plt
from src.pack import Pack

def plot(df, pack: Pack):

    with plt.style.context("dark_background"):
        fig, ax1 = plt.subplots()

        ax1.set_xlabel('Time (s)')
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        colors = set(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        temps = colors.copy()

        for j in range(pack.parallel):
            if (len(temps) == 0):
                temps = colors.copy()

            c = temps.pop()

            ax1.set_ylabel('Pack Voltage (V)', color='tab:red')
            voltname = pack.volts[j].name # cells[i,j].volt.name
            ax1.plot(df['Time'], df[voltname], color=c, label=voltname)
            ax1.tick_params(axis='y', labelcolor='tab:red')

            ax2.set_ylabel('String Currents (A/m2)', color='tab:blue')  # we already handled the x-label with ax1
            iname = pack.iapps[j].name
            ax2.plot(df['Time'], df[iname], color=c, label=iname)
            ax2.tick_params(axis='y', labelcolor='tab:blue')


        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0 + box.height * 0.2,
                        box.width, box.height * 0.8])

        # Put a legend below current axis
        ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
                fancybox=True, shadow=True, ncol=2)

        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0 + box.height * 0.2,
                        box.width, box.height * 0.8])

        # Put a legend below current axis
        ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.45),
                fancybox=True, shadow=True, ncol=2)

        plt.grid(linewidth=0.2)
        plt.show()