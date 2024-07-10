import pandas as pd
from matplotlib import pyplot as plt
import params as p

def plotter(celldf: pd.DataFrame):
    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Time (s)')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax1.set_ylabel('Length (m)', color='black')
    ax1.plot(celldf['Time'], celldf['Cell1 Anode vSEI Length'], color='black')
    ax1.tick_params(axis='y', labelcolor='black')

    # ax1.set_ylabel('Concentration', color='green')
    # ax1.plot(celldf['Time'], celldf['Cell1 Cathode vConcentration'], color='green')
    # ax1.plot(celldf['Time'], celldf['Cell1 Anode vConcentration'], color='green')
    # ax1.tick_params(axis='y', labelcolor='green')

    # ax2.set_ylabel('Voltage (V)', color='tab:blue')  # we already handled the x-label with ax1
    # ax2.plot(celldf['Time'], celldf['Cell1 Cathode vPhi'] - celldf['Cell1 Anode vPhi'], color='blue')
    # ax2.tick_params(axis='y', labelcolor='tab:blue')

    plt.grid(linewidth=0.2)
    plt.tight_layout()
    plt.legend()
    plt.show()