import numpy as np
import pandas as pd 
from matplotlib import pyplot as plt
import json
from src.pack import Pack
from src.cell import Cell
import sys
import os
import pickle

DISCHARGE = "CC-discharge"
CV_CHARGE = "CV-charge"
CC_CHARGE = "CC-charge"
CHARGE = "CV-charge|CC-charge"

class Experiment:
    def __init__(self, experiment: str):
        self.path = f"data/{experiment}/"

        sys.path.append(os.path.join(os.getcwd(), "src"))

        self.pack = None
        with open(self.path+"model.pkl", 'rb') as f:
            self.pack = pickle.load(f)

        self.profile = None
        with open(self.path+"profile.json", 'r') as f:
            self.profile = json.load(f)

        self.profile_str =json.dumps(self.profile, indent=4)

        self.data = pd.read_csv(self.path+"data.csv", index_col=[0,1,2])
        self.caps = pd.read_csv(self.path+"capacities.csv", index_col=0)
        self.CYCLE = self.data.index.get_level_values(0)
        self.PROTOCOL = self.data.index.get_level_values(1)

    def __str__(self):
        return self.profile_str

    def select_cycles(self, cycles=[], protocols=[]):
        CYCLE = self.data.index.get_level_values(0)
        PROTOCOL = self.data.index.get_level_values(1)

        flts = [True, True]
        if len(cycles) != 0:
            flts[0] = CYCLE.isin(cycles)
            #self.data = self.data.loc[self.CYCLE.isin(cycles)]

        if len(protocols) != 0:
            flts[1] = PROTOCOL.str.contains("|".join(protocols))
            # print("|".join(protocols))
            # self.data = self.data.loc[
                # self.PROTOCOL.str.contains("|".join(protocols))
            # ]

        self.data = self.data.loc[flts[0] & flts[1]]

    def select_attributes(self, attrs: list):
        joined = '|'.join(attrs)
        self.data = self.data.filter(regex=f'Time|{joined}')

    def plotter(self, isolate_cycles=True):
        # Helper function to encapsulate the plotting logic
        def plot_columns(data, t, label_prefix=''):
            """Helper function to plot columns."""
            for col in data.columns.drop(['Time', 'Global Time']):
                label = f'{label_prefix}{col}' if label_prefix else col
                plt.plot(data[t], data[col], label=label)
        
        if isolate_cycles:
            t = 'Time'
            for cnum, group in self.data.groupby(level=0):
                plot_columns(group, t, label_prefix=f'C{cnum}_')
        else:
            t = 'Global Time'
            plot_columns(self.data, t)

        plt.xlabel(t)
        plt.legend()
        plt.show()

    def plot_capacities(self, cycles=[], cells=[]):
        cyc = self.caps.index
        cap_data = self.caps
        if len(cycles) != 0:
            cyc = cycles
            cap_data = self.caps.loc[cycles]

        cell_list = self.caps.columns
        if type(cells) is np.ndarray:
            # 'cells' looks like self.pack.cells[IND: IND]
            cell_list = [cell.name for cell in cells.flatten()]

        elif type(cells) is Cell:
            cell_list = [cells.name]

        fig, ax = plt.subplots()

        # Plot data with a legend
        for column in cell_list:
            ax.scatter(cyc, cap_data[column], label=column)

        # Customize labels and title
        ax.set_xlabel('Cycle #')
        ax.set_ylabel('Discharge Capacity (A/m2)')
        ax.set_title('Discharge Capacity by Cycle #')

        # Move the legend (just annoying configuration)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.2,
                        box.width, box.height * 0.8])

        ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4), 
                ncol=self.pack.series, fancybox=True, shadow=True)

        # Add grid and show plot
        ax.grid()
        plt.show()

    def get_capacities(self) -> pd.DataFrame:
        return self.caps

    def get_data(self) -> pd.DataFrame:
        return self.data

    def get_pack(self) -> Pack:
        return self.pack

    def get_profile(self) -> dict:
        return self.profile

    def reset(self) -> None:
        self.data = pd.read_csv(self.path+"data.csv", index_col=[0,1,2])

    def to_csv(self, filename: str) -> None:
        self.data.to_csv(self.path+filename, index=True)
