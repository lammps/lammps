from abc import ABC, abstractmethod
import pickle

class MLIAPUnified(ABC):
    """Abstract base class for MLIAPUnified."""

    def __init__(self):
        self.interface = None
        self.element_types = None
        self.ndescriptors = None
        self.nparams = None
        self.rcutfac = None

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
