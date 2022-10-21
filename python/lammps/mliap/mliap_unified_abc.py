from abc import ABC, abstractmethod
import pickle

class MLIAPUnified(ABC):
    """Abstract base class for MLIAPUnified."""

    def __init__(self, interface = None, element_types = None,
                 ndescriptors = None, nparams = None, rcutfac = None):
        self.interface = interface
        self.element_types = element_types
        self.ndescriptors = ndescriptors
        self.nparams = nparams
        self.rcutfac = rcutfac

    @abstractmethod
    def compute_gradients(self, data):
        """Compute gradients."""

    @abstractmethod
    def compute_descriptors(self, data):
        """Compute descriptors."""

    @abstractmethod
    def compute_forces(self, data):
        """Compute forces."""

    def pickle(self, fname):
        with open(fname, 'wb') as fp:
            pickle.dump(self, fp)
