from multiprocessing.sharedctypes import Value
import sys
import pandas as pd 
from matplotlib import pyplot as plt
import numpy as np
import pickle
import json
from src.pack import Pack

DISCHARGE = "CC-discharge"
CV_CHARGE = "CV-charge"
CC_CHARGE = "CC-charge"
CHARGE = "CV-charge|CC-charge"

class Experiment:
    def __init__(self, experiment: str):
        self.path = f"data/{experiment}/"

        profile = None
        with open(self.path+"profile.json", 'r') as f:
            profile = json.load(f)

        self.profile =json.dumps(profile, indent=4)

        self.data = pd.read_csv(self.path+"data.csv", index_col=[0,1,2])
        self.CYCLE = self.data.index.get_level_values(0)
        self.PROTOCOL = self.data.index.get_level_values(1)

    def show_profile(self):
        print(self.profile)

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

    def plotter(self, columns=[], local_time=False):
        t = 'Global Time'
        if local_time:
            t = 'Time'
    
        d = self.data
        if len(columns) != 0:
            pass
            # d = self.data.loc[: columns]
            # joined = '|'.join(columns)
            # a = d.filter(regex=f"{t}|{joined}")
            # return lambda: a.plot(x=t, y=d.filter(regex=joined).columns, kind='line')

        else:
            for col in d.columns.drop(['Time', 'Global Time']):
                plt.plot(d[t], d[col], label=col)

            plt.xlabel(t)
            plt.legend()
            plt.show()
            # return lambda: d.plot(x=t, y=d.columns.drop(t), kind='line')

    def get_data(self):
        return self.data

    def plot(self):
        # Plot the experiment data
        pass

    def reset(self):
        self.data = pd.read_csv(self.path+"data.csv", index_col=[0,1,2])
