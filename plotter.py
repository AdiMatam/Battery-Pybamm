from tkinter import W
from  matplotlib import pyplot as plt
from typing import List
from cell import Cell

def plot(df, cells: List[Cell]):
    # seaborn-talk
    with plt.style.context("seaborn-v0_8-talk"):
        fig, ax1 = plt.subplots()

        ax1.set_xlabel('Time (s)')
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        colors = set(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        temps = colors.copy()

        for i in range(len(cells)):
            cell = cells[i]
            if (len(temps) == 0):
                temps = colors.copy()

            c = temps.pop()

            ax1.set_ylabel('Concentration (mol / m2)', color='tab:red')
            ax1.plot(df['Time'], df[cell.pos.surf_csn_name], color=c, label=cell.pos.surf_csn_name)
            ax1.plot(df['Time'], df[cell.neg.surf_csn_name], color=c, label=cell.neg.surf_csn_name)
            ax1.tick_params(axis='y', labelcolor='tab:red')

            ax2.set_ylabel('Potential (V)', color='tab:blue')  # we already handled the x-label with ax1
            ax2.plot(df['Time'], df[cell.voltage_name], color=c, label=cell.voltage_name)
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

        plt.show()