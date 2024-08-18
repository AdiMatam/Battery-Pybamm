from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE

smallpack = Experiment("1by2_5cycles_isocz")
#smallpack.show_profile()

smallpack.select_cycles(
    cycles=[1, 5], 
    protocols=[CC_CHARGE]
)

smallpack.select_attributes(["SEI"])

print(smallpack.get_data())
# d = smallpack.get_data()
# for key, group in d.groupby(level=0):
    # print(f"Section: {key}")
    # print(group)
    # print()

smallpack.plotter(local_time=True)

# smallpack.reset()
# print(smallpack.get_data())