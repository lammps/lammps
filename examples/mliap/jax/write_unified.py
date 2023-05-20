"""
interface for creating LAMMPS MLIAP Unified models.
"""
import pickle

import numpy as np

from lammps.mliap.mliap_unified_abc import MLIAPUnified
#from deploy_script import MyModel

class MLIAPInterface(MLIAPUnified):
    """
    Class for creating ML-IAP Unified model based on hippynn graphs.
    """
    def __init__(self, model, element_types, cutoff=4.5, ndescriptors=1):
        """
        :param model: class defining the model
        :param element_types: list of atomic symbols corresponding to element types
        :param ndescriptors: the number of descriptors to report to LAMMPS
        :param model_device: the device to send torch data to (cpu or cuda)
        """
        super().__init__()
        self.model = model
        self.element_types = element_types
        self.ndescriptors = ndescriptors
        #self.model_device = model_device
        

        # Build the calculator
        # TODO: Make this cutoff depend on model cutoff, ideally from deployed model itself but could 
        # be part of deploy step.
        #rc = 4.5
        self.rcutfac = 0.5*cutoff # Actual cutoff will be 2*rc
        #print(self.model.nparams)
        self.nparams = 10
        #self.rcutfac, self.species_set, self.graph = setup_LAMMPS()
        #self.nparams = sum(p.nelement() for p in self.graph.parameters())
        #self.graph.to(torch.float64)

    def compute_descriptors(self, data):
        pass

    def compute_gradients(self, data):
        pass

    def compute_forces(self, data):
        #print(">>>>> hey!")
        #elems = self.as_tensor(data.elems).type(torch.int64).reshape(1, data.ntotal)

        """
        elems = self.as_tensor(data.elems).type(torch.int64) + 1
        #z_vals = self.species_set[elems+1]
        pair_i = self.as_tensor(data.pair_i).type(torch.int64)
        pair_j = self.as_tensor(data.pair_j).type(torch.int64)
        rij = self.as_tensor(data.rij).type(torch.float64).requires_grad_(True)
        nlocal = self.as_tensor(data.nlistatoms) 
        """

        rij = data.rij

        #(total_energy, fij) = self.network(rij, None, None, None, nlocal, elems, pair_i, pair_j, "cpu", dtype=torch.float64, mode="lammps")

        test = self.model(rij)
         
        #data.update_pair_forces(fij)
        #data.energy = total_energy.item()
 
        pass

def setup_LAMMPS(energy):
    """

    :param energy: energy node for lammps interface
    :return: graph for computing from lammps MLIAP unified inputs.
    """

    model = TheModelClass(*args, **kwargs)

    save_state_dict = torch.load("Ta_Pytorch.pt")
    model.load_state_dict(save_state_dict["model_state_dict"])


    #model.load_state_dict(torch.load(PATH))
    model.eval()
    
    #model.eval()
    return model