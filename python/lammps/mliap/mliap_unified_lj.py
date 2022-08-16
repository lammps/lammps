from .mliap_unified_abc import MLIAPUnified
import numpy as np


class MLIAPUnifiedLJ(MLIAPUnified):
    """Test implementation for MLIAPUnified."""

    def __init__(self):
        super().__init__()
        self.element_types = ["Ar"]
        self.ndescriptors = 1
        self.nparams = 3
        # Mimicking the LJ pair-style:
        # pair_style lj/cut 2.5
        # pair_coeff * * 1 1
        self.epsilon = 1.0
        self.sigma = 1.0
        self.rcutfac = 1.25

    def compute_gradients(self, data):
        """Test compute_gradients."""
    
    def compute_descriptors(self, data):
        """Test compute_descriptors."""
    
    def compute_forces(self, data):
        """Test compute_forces."""
        eij, fij = self.compute_pair_ef(data)
        data.update_pair_energy(eij)
        data.update_pair_forces(fij)
 
    def compute_pair_ef(self, data):
        rij = data.rij

        r2inv = 1.0 / np.sum(rij ** 2, axis=1)
        r6inv = r2inv * r2inv * r2inv

        lj1 = 4.0 * self.epsilon * self.sigma**12
        lj2 = 4.0 * self.epsilon * self.sigma**6

        eij = r6inv * (lj1 * r6inv - lj2)
        fij = r6inv * (3.0 * lj2 - 6.0 * lj2 * r6inv) * r2inv
        fij = fij[:, np.newaxis] * rij
        return eij, fij
