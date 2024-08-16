from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE

smallpack = Experiment("1by2_5cycles_isocz")
#smallpack.show_profile()

smallpack.select_cycles(
    cycles=[1, 5], 
    protocols=[CC_CHARGE]
)

smallpack.select_attributes(["Cathode Concentration"])

#print(smallpack.get_data())
smallpack.plotter(local_time=True)

# smallpack.reset()
# print(smallpack.get_data())