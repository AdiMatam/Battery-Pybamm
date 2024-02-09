from  matplotlib import pyplot as plt

def plot(df, cells):
    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Time (s)')

    color = 'tab:red'
    ax1.set_ylabel('Concentration (mol / m2)', color=color)
    ax1.plot(df['Time'], df[cells[0].pos.surf_csn_name], color=color)
    ax1.plot(df['Time'], df[cells[0].neg.surf_csn_name], color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('Potential (V)', color=color)  # we already handled the x-label with ax1
    ax2.plot(df['Time'], df[cells[0].voltage_name], color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    plt.tight_layout()
    plt.show()