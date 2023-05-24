from lammps.mliap.mliap_unified_abc import MLIAPUnified
import numpy as np
import jax
import jax.numpy as jnp
from jax import jit
from functools import partial


class MLIAPUnifiedJAX(MLIAPUnified):
    """Test implementation for MLIAPUnified."""

    def __init__(self, element_types, epsilon=1.0, sigma=1.0, rcutfac=1.25):
        # ARGS: interface, element_types, ndescriptors, nparams, rcutfac
        super().__init__(None, element_types, 1, 3, rcutfac)
        # Mimicking the LJ pair-style:
        # pair_style lj/cut 2.5
        # pair_coeff * * 1 1
        self.epsilon = epsilon
        self.sigma = sigma
        self.npair_max = 250000

    def compute_gradients(self, data):
        """Test compute_gradients."""

    def compute_descriptors(self, data):
        """Test compute_descriptors."""

    def compute_forces(self, data):
        """Test compute_forces."""
        rij  = data.rij

        # TODO: Take max npairs from the LAMMPS Cython side.
        if (data.npairs > self.npair_max):
            self.npair_max = data.npairs

        npad = self.npair_max - data.npairs
        # TODO: Take pre-padded rij from the LAMMPS Cython side.
        #       This might account for ~2-3x slowdown compared to original LJ.
        rij = np.pad(rij, ((0,npad), (0,0)), 'constant')

        eij, fij = self.compute_pair_ef(rij)

        data.update_pair_energy(np.array(np.double(eij)))
        data.update_pair_forces(np.array(np.double(fij)))

    #@jax.jit # <-- This will error! See https://github.com/google/jax/issues/1251
    # @partial takes a function (e.g. jax.jit) as an arg.
    @partial(jax.jit, static_argnums=(0,))
    def compute_pair_ef(self, rij):

        r2inv = 1.0 / np.sum(rij ** 2, axis=1)
        r6inv = r2inv * r2inv * r2inv

        lj1 = 4.0 * self.epsilon * self.sigma**12
        lj2 = 4.0 * self.epsilon * self.sigma**6

        eij = r6inv * (lj1 * r6inv - lj2)
        fij = r6inv * (3.0 * lj2 - 6.0 * lj2 * r6inv) * r2inv
        fij = fij[:, np.newaxis] * rij
        return eij, fij