import numpy as np
import pickle
from pathlib import Path
from write_unified import MLIAPInterface

class MyModel():
  def __init__(self,blah):
      """
      coeffs = np.genfromtxt(file,skip_header=6)
      self.bias = coeffs[0]
      self.weights = coeffs[1:]
      """
      self.blah = blah
      self.n_params = 3 #len(coeffs)
      self.n_descriptors = 1 #len(self.weights)
      self.n_elements = 1

  def __call__(self,rij):
      print(rij)
      #energy[:] = bispectrum @ self.weights + self.bias
      #beta[:] = self.weights
      return 5

model = MyModel(1)

#unified = MLIAPInterface(model, ["Ta"], model_device="cpu")

def create_pickle():
    unified = MLIAPInterface(model, ["Ta"])
    unified.pickle('mliap_jax.pkl')

create_pickle()