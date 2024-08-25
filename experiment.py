import pandas as pd 
from matplotlib import pyplot as plt
import json
from src.pack import Pack

DISCHARGE = "CC-discharge"
CV_CHARGE = "CV-charge"
CC_CHARGE = "CC-charge"
CHARGE = "CV-charge|CC-charge"

class Experiment:
    def __init__(self, experiment: str):
        self.path = f"data/{experiment}/"

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
        flts = [True, True]
        if len(cycles) != 0:
            flts[0] = self.CYCLE.isin(cycles)
            #self.data = self.data.loc[self.CYCLE.isin(cycles)]

        if len(protocols) != 0:
            flts[1] = self.PROTOCOL.str.contains("|".join(protocols))
            # print("|".join(protocols))
            # self.data = self.data.loc[
                # self.PROTOCOL.str.contains("|".join(protocols))
            # ]

        self.data = self.data.loc[flts[0] & flts[1]]

    def select_attributes(self, attrs: list):
        joined = '|'.join(attrs)
        self.data = self.data.filter(regex=f'Time|{joined}')

    def plotter(self, local_time=False, isolate_cycles=True):
        # Helper function to encapsulate the plotting logic
        def plot_columns(data, t, label_prefix=''):
            """Helper function to plot columns."""
            for col in data.columns.drop(['Time', 'Global Time']):
                label = f'{label_prefix}{col}' if label_prefix else col
                plt.plot(data[t], data[col], label=label)
        

        t = 'Global Time'
        if local_time:
            t = 'Time'

        if isolate_cycles:
            for cnum, group in self.data.groupby(level=0):
                plot_columns(group, t, label_prefix=f'C{cnum}_')
        else:
            plot_columns(self.data, t)

        plt.xlabel(t)
        plt.legend()
        plt.show()

    def plot_capacities(self):
        for column in self.caps.columns:
            plt.scatter(self.caps.index, self.caps[column], label=column)

        plt.xlabel('Cycle #')
        plt.ylabel('Discharge Capacity (A/m2)')
        plt.title('Discharge Capacity by Cycle #')

        plt.legend()
        plt.show()
        self.caps

    def get_capacities(self):
        return self.caps

    def get_data(self):
        return self.data

    def reset(self):
        self.data = pd.read_csv(self.path+"data.csv", index_col=[0,1,2])

    def to_csv(self, filename: str):
        self.data.to_csv(self.path+filename, index=True)
