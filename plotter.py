from tkinter import W
from  matplotlib import pyplot as plt
import random
from matplotlib import colors as mcolors
from typing import List
from cell import Cell

def generate_random_color():
    while True:
        r, g, b = random.random(), random.random(), random.random()
        color = "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))
        # Check if the color is too light for a white background
        if mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])[2] < 0.4:
            return color


def plot(df, cells: List[Cell]):
    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Time (s)')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    for cell in cells:
        color = 'tab:red'
        c = generate_random_color()
        ax1.set_ylabel('Concentration (mol / m2)', color=color)
        ax1.plot(df['Time'], df[cell.pos.surf_csn_name], color=c, label=cell.pos.surf_csn_name)
        ax1.plot(df['Time'], df[cell.neg.surf_csn_name], color=c, label=cell.neg.surf_csn_name)
        ax1.tick_params(axis='y', labelcolor=color)

        color = 'tab:blue'
        ax2.set_ylabel('Potential (V)', color=color)  # we already handled the x-label with ax1
        ax2.plot(df['Time'], df[cell.voltage_name], color=generate_random_color(), label=cell.voltage_name)
        ax2.tick_params(axis='y', labelcolor=color)

    plt.legend()
    plt.tight_layout()
    plt.show()