# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   https://www.lammps.org/ Sandia National Laboratories
#   LAMMPS Development team: developers@lammps.org
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

# ----------------------------------------------------------------------
#   Contributing author: Nicholas Lubbers (LANL)
# -------------------------------------------------------------------------

import numpy as np
import torch

def calc_n_params(model):
    """
    Returns the sum of two decimal numbers in binary digits.

        Parameters:
                model (torch.nn.Module): Network model that maps descriptors to a per atom attribute

        Returns:
                n_params (int): Number of NN model parameters
    """
    return sum(p.nelement() for p in model.parameters())

class TorchWrapper(torch.nn.Module):
    """
    A class to wrap Modules to ensure lammps mliap compatability.

    ...

    Attributes
    ----------
    model : torch.nn.Module
        Network model that maps descriptors to a per atom attribute

    device : torch.nn.Module (None)
        Accelerator device

    dtype : torch.dtype (torch.float64)
        Dtype to use on device

    n_params : torch.nn.Module (None)
        Number of NN model parameters

    n_descriptors : int
        Max number of per atom descriptors

    n_elements : int
        Max number of elements


    Methods
    -------
    forward(descriptors, elems):
        Feeds descriptors to network model to produce per atom energies and forces.
    """

    def __init__(self, model, n_descriptors, n_elements, n_params=None, device=None, dtype=torch.float64):
        """
        Constructs all the necessary attributes for the network module.

        Parameters
        ----------
            model : torch.nn.Module
                Network model that maps descriptors to a per atom attribute

            n_descriptors : int
                Max number of per atom descriptors

            n_elements : int
                Max number of elements

            n_params : torch.nn.Module (None)
                Number of NN model parameters

            device : torch.nn.Module (None)
                Accelerator device

            dtype : torch.dtype (torch.float64)
                Dtype to use on device
        """

        super().__init__()
        self.model = model
        self.device = device
        self.dtype = dtype

        # Put model on device and convert to dtype
        self.to(self.dtype)
        self.to(self.device)

        if n_params is None:
            n_params = calc_n_params(model)

        self.n_params = n_params
        self.n_descriptors = n_descriptors
        self.n_elements = n_elements

    def forward(self, elems, descriptors, beta, energy,use_gpu_data=False):
        """
        Takes element types and descriptors calculated via lammps and
        calculates the per atom energies and forces.

        Parameters
        ----------
        elems : numpy.array
            Per atom element types

        descriptors : numpy.array
            Per atom descriptors

        beta : numpy.array
            Expired beta array to be filled with new betas

        energy : numpy.array
            Expired per atom energy array to be filled with new per atom energy
            (Note: This is a pointer to the lammps per atom energies)


        Returns
        -------
        None
        """
        descriptors = torch.as_tensor(descriptors,dtype=self.dtype, device=self.device).requires_grad_(True)
        elems = torch.as_tensor(elems,dtype=torch.int32, device=self.device)
        elems=elems-1
        with torch.autograd.enable_grad():

            if (use_gpu_data):
                energy_nn = torch.as_tensor(energy,dtype=self.dtype, device=self.device)
                energy_nn[:] = self.model(descriptors, elems).flatten()
            else:
                energy_nn = self.model(descriptors, elems).flatten()
                energy[:] = energy_nn.detach().cpu().numpy().astype(np.float64)

        if (use_gpu_data):
            beta_nn = torch.as_tensor(beta,dtype=self.dtype, device=self.device)
            beta_nn[:] = torch.autograd.grad(energy_nn.sum(), descriptors)[0]
        else:
            beta_nn = torch.autograd.grad(energy_nn.sum(), descriptors)[0]
            beta[:] = beta_nn.detach().cpu().numpy().astype(np.float64)

class IgnoreElems(torch.nn.Module):
    """
    A class to represent a NN model agnostic of element typing.

    ...

    Attributes
    ----------
    subnet : torch.nn.Module
        Network model that maps descriptors to a per atom attribute

    Methods
    -------
    forward(descriptors, elems):
        Feeds descriptors to network model
    """

    def __init__(self, subnet):
        """
        Constructs all the necessary attributes for the network module.

        Parameters
        ----------
            subnet : torch.nn.Module
                Network model that maps descriptors to a per atom attribute
        """

        super().__init__()
        self.subnet = subnet

    def forward(self, descriptors, elems):
        """
        Feeds descriptors to network model

        Parameters
        ----------
        descriptors : torch.tensor
            Per atom descriptors

        elems : torch.tensor
            Per atom element types

        Returns
        -------
        self.subnet(descriptors) : torch.tensor
            Per atom attribute computed by the network model
        """

        return self.subnet(descriptors)


class UnpackElems(torch.nn.Module):
    """
    A class to represent a NN model pseudo-agnostic of element typing for
    systems with multiple element typings.

    ...

    Attributes
    ----------
    subnet : torch.nn.Module
        Network model that maps descriptors to a per atom attribute

    n_types : int
        Number of atom types used in training the NN model.

    Methods
    -------
    forward(descriptors, elems):
        Feeds descriptors to network model after adding zeros into
        descriptor columns relating to different atom types
    """

    def __init__(self, subnet, n_types):
        """
        Constructs all the necessary attributes for the network module.

        Parameters
        ----------
            subnet : torch.nn.Module
                Network model that maps descriptors to a per atom attribute.

            n_types : int
                Number of atom types used in training the NN model.
        """
        super().__init__()
        self.subnet = subnet
        self.n_types = n_types

    def forward(self, descriptors, elems):
        """
        Feeds descriptors to network model after adding zeros into
        descriptor columns relating to different atom types

        Parameters
        ----------
        descriptors : torch.tensor
            Per atom descriptors

        elems : torch.tensor
            Per atom element types

        Returns
        -------
        self.subnet(descriptors) : torch.tensor
            Per atom attribute computed by the network model
        """

        unpacked_descriptors = torch.zeros(elems.shape[0], self.n_types, descriptors.shape[1], dtype=torch.float64)
        for i, ind in enumerate(elems):
            unpacked_descriptors[i, ind, :] = descriptors[i]
        return self.subnet(torch.reshape(unpacked_descriptors, (elems.shape[0], -1)), elems)


class ElemwiseModels(torch.nn.Module):
    """
    A class to represent a NN model dependent on element typing.

    ...

    Attributes
    ----------
    subnets : list of torch.nn.Modules
        Per element type network models that maps per element type
        descriptors to a per atom attribute.

    n_types : int
        Number of atom types used in training the NN model.

    Methods
    -------
    forward(descriptors, elems):
        Feeds descriptors to network model after adding zeros into
        descriptor columns relating to different atom types
    """

    def __init__(self, subnets, n_types):
        """
        Constructs all the necessary attributes for the network module.

        Parameters
        ----------
            subnets : list of torch.nn.Modules
                Per element type network models that maps per element
                type descriptors to a per atom attribute.

            n_types : int
                Number of atom types used in training the NN model.
        """

        super().__init__()
        self.subnets = subnets
        self.n_types = n_types

    def forward(self, descriptors, elems, dtype=torch.float64):
        """
        Feeds descriptors to network model after adding zeros into
        descriptor columns relating to different atom types

        Parameters
        ----------
        descriptors : torch.tensor
            Per atom descriptors

        elems : torch.tensor
            Per atom element types

        Returns
        -------
        self.subnets(descriptors) : torch.tensor
            Per atom attribute computed by the network model
        """

        self.dtype=dtype
        self.to(self.dtype)

        per_atom_attributes = torch.zeros(elems.size(dim=0), dtype=self.dtype)
        given_elems, elem_indices = torch.unique(elems, return_inverse=True)
        for i, elem in enumerate(given_elems):
            self.subnets[elem].to(self.dtype)
            per_atom_attributes[elem_indices == i] = self.subnets[elem](descriptors[elem_indices == i]).flatten()
        return per_atom_attributes
