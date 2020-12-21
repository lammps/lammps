import sys
import numpy as np
import torch

# torch.nn.modules useful for defining a MLIAPPY model.
from lammps.mliap.pytorch import TorchWrapper, IgnoreElems

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

# Wrap the pytorch model for usage with mliappy coupling.
model = IgnoreElems(lin) # The linear module does not use the types.
n_descriptors = lin.weight.shape[1]
n_elements = 1
linked_model = TorchWrapper(model,n_descriptors=n_descriptors,n_elements=n_elements)

torch.save(linked_model,"Ta06A.mliap.pytorch.model.pt")
