import pickle

with open("cells.pkl", "rb") as f:
    packarray = pickle.load(f)

    print(packarray.shape)
    print(packarray[0,0].pos.c.name)