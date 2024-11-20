import pandas as pd
from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE

## LOADING AN EXPERIMENT (folder name within data/)
squarepack = Experiment("5by5_1.0C")
PACK = squarepack.get_pack()

## PRINTS THE `PROFILE.JSON` as string (operating condition data)
# print(squarepack)

## DATA FILTERING, SELECTION

### Selection of cycles and/or protocols. Both lists are OPTIONAL 
### Corresponds to rows in the experiment/data.csv
squarepack.select_cycles(
    cycles=[2], 
    protocols=[DISCHARGE]
)

### Selection of cell attributes (columns)
### Uses 'fuzzy' regex searching -- so searching for ['SEI'] will choose ALL columns with SEI in the name
### In other words, each cell's SEI data will be in the filtered table
# squarepack.select_attributes(["SEI"])
squarepack.select_attributes(["Cell.*?Voltage"])

df = squarepack.get_data().copy(deep=True)
voltages = df.filter(regex="Cell.*?Voltage")
lower_bound = PACK.voltage_window[0] / PACK.series
last_volts = voltages.iloc[-1]
filtered = last_volts[last_volts < lower_bound]

result_df = pd.DataFrame({
    'Column': filtered.index,
    'OD-Voltage': filtered.values
})

print(result_df)


### Plot the CURRENT dataset (i.e. after all predecessing filters)
### isolate_cycles =True:  Plot data for EACH cycle as separate line with respect to "local time"
###                =False: Plot data with respect to "global time" (no delineation by cycle #)
squarepack.plotter(isolate_cycles=True)

# ### To get the dataframe. Modifications can be made outside for those well-versed in pandas...
# df = squarepack.get_data()

# ### The various functions will manipulate the 'dataframe'. In order to reset to the contents of the file (remove all filterations)
# squarepack.reset()


# ## (MIGHT WORK? STILL PENDING DEBUG) READING DISCHARGE CAPACITY DATA
# squarepack.plot_capacities(cycles=[1] + list(range(10, PACK.cycles, 20)), strings=[0,1])