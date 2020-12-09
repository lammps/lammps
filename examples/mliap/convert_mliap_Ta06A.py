import sys
import numpy as np
import torch
import pickle
import os
import shutil

shutil.copyfile('../../src/MLIAP/mliappy_pytorch.py','./mliappy_pytorch.py')

import mliappy_pytorch

# Read coefficients
coeffs = np.genfromtxt("Ta06A.mliap.model",skip_header=6)

# Write coefficients to a pytorch linear model
bias = coeffs[0]
weights = coeffs[1:]
lin = torch.nn.Linear(weights.shape[0],1)
lin.to(torch.float64)
with torch.autograd.no_grad():
    lin.weight.set_(torch.from_numpy(weights).unsqueeze(0))
    lin.bias.set_(torch.as_tensor(bias,dtype=torch.float64).unsqueeze(0))

# Wrap the pytorch model for usage with mliappy energy model
model = mliappy_pytorch.IgnoreElems(lin)
n_descriptors = lin.weight.shape[1]
n_params = mliappy_pytorch.calc_n_params(model)
n_types = 1
linked_model = mliappy_pytorch.TorchWrapper64(model,n_descriptors=n_descriptors,n_elements=n_types)

# Save the result.
with open("Ta06A.mliap.pytorch.model.pkl",'wb') as pfile:
    pickle.dump(linked_model,pfile)
