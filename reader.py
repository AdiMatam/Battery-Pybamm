from experiment import Experiment, CHARGE, CC_CHARGE, CV_CHARGE, DISCHARGE

squarepack = Experiment("5by5_100cycles_const")
print(squarepack)

squarepack.plot_capacities()

# squarepack.select_cycles(
    # cycles=[1, 5], 
    # protocols=[CC_CHARGE]
# )

# squarepack.select_attributes(["SEI"])

# print(squarepack.get_data())
# # d = smallpack.get_data()
# # for key, group in d.groupby(level=0):
    # # print(f"Section: {key}")
    # # print(group)
    # # print()

# squarepack.plotter(local_time=True)

# # smallpack.reset()
# # print(smallpack.get_data())