import pandas as pd
from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE
from matplotlib import pyplot as plt

## LOADING AN EXPERIMENT (folder name within data/)
squarepack = Experiment("Single_1.0C_3.0")
PACK = squarepack.get_pack()

## PRINTS THE `PROFILE.JSON` as string (operating condition data)
# print(squarepack)

## DATA FILTERING, SELECTION

### Selection of cycles and/or protocols. Both lists are OPTIONAL 
### Corresponds to rows in the experiment/data.csv
squarepack.select_cycles(
    cycles=[1,11,21,31,41], #+ list(range(10,PACK.cycles, 10)),
    protocols=[CC_CHARGE, CV_CHARGE]
)

### Selection of cell attributes (columns)
### Uses 'fuzzy' regex searching -- so searching for ['SEI'] will choose ALL columns with SEI in the name
### In other words, each cell's SEI data will be in the filtered table
# squarepack.select_attributes(["SEI"])
squarepack.select_attributes(["Side"])

### Plot the CURRENT dataset (i.e. after all predecessing filters)
### isolate_cycles =True:  Plot data for EACH sub-cycle as separate line with respect to "local time"
###                =False: Plot data with respect to "global time" (no delineation by cycle #)

squarepack.plotter(isolate_cycles=True)

df = squarepack.get_data().copy(deep=True)
print(len(df.iloc[:,-1].tolist()))
print(*df.iloc[:,-1].tolist(), sep='\n', end='\n\n')

plt.legend()
plt.show()

