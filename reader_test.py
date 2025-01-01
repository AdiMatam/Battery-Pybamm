import pandas as pd
from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE
from matplotlib import pyplot as plt

## LOADING AN EXPERIMENT (folder name within data/)
squarepack = Experiment("Single_1.0C_3.0_noage")
PACK = squarepack.get_pack()

## PRINTS THE `PROFILE.JSON` as string (operating condition data)
# print(squarepack)

## DATA FILTERING, SELECTION

### Selection of cycles and/or protocols. Both lists are OPTIONAL 
### Corresponds to rows in the experiment/data.csv
squarepack.select_cycles(
    cycles=[2], #+ list(range(10,PACK.cycles, 10)),
)

### Selection of cell attributes (columns)
### Uses 'fuzzy' regex searching -- so searching for ['SEI'] will choose ALL columns with SEI in the name
### In other words, each cell's SEI data will be in the filtered table
# squarepack.select_attributes(["SEI"])
squarepack.select_attributes(["Cathode.*?Concentration"])

### OVER-DISCHARGE CHECK
# df = squarepack.get_data().copy(deep=True)
# voltages = df.filter(regex="Cell.*?Voltage")
# lower_bound = PACK.voltage_window[0] / PACK.series
# last_volts = voltages.iloc[-1]
# filtered = last_volts[last_volts < lower_bound]

# result_df = pd.DataFrame({
#     'Column': filtered.index,
#     'OD-Voltage': filtered.values
# })

# print(result_df)

### Plot the CURRENT dataset (i.e. after all predecessing filters)
### isolate_cycles =True:  Plot data for EACH sub-cycle as separate line with respect to "local time"
###                =False: Plot data with respect to "global time" (no delineation by cycle #)

from matlab_export import t,cathode

squarepack.plotter(isolate_cycles=True)
plt.plot(t, cathode, label="Matlab Code")

df = squarepack.get_data().copy(deep=True)
cat_python = df.iloc[:,-1].tolist()

#print(df['Time'][-1] - t[-1])
print(cat_python[0] - cathode[0])
print(cat_python[-1] - cathode[-1])


plt.legend()
plt.show()
