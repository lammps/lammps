import pickle
from mytest import MLIAPMod

mymodel = MLIAPMod(["H", "C", "O"])
with open("my_model.pkl", "wb") as f:
    pickle.dump(mymodel, f)
