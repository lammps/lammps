from lammps.mliap.mliap_unified_abc import MLIAPUnified
import numpy as np
import jax
import jax.dlpack
import jax.numpy as jnp
from jax import jit
from functools import partial
import cupy
import os

# Required else get `jaxlib.xla_extension.XlaRuntimeError: RESOURCE_EXHAUSTED: Out of memory`
# Does not fix GPU problem with larger num. atoms.
#os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"]="false"
#os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"]=".XX"
#os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"]="platform"

@jax.jit
def lj_potential(epsilon, sigma, rij):
    # A pure function we can differentiate:
    def _tot_e(rij):
        r2inv = 1.0 / jnp.sum(rij ** 2, axis=1)
        r6inv = r2inv * r2inv * r2inv

        lj1 = 4.0 * epsilon * sigma**12
        lj2 = 4.0 * epsilon * sigma**6

        eij = r6inv * (lj1 * r6inv - lj2)
        return 0.5 * jnp.sum(eij), eij
    # Construct a function computing _tot_e and its derivative
    (_, eij), fij = jax.value_and_grad(_tot_e, has_aux=True)(rij)
    return eij, fij


class MLIAPUnifiedJAXKokkos(MLIAPUnified):
    """JAX wrapper for MLIAPUnified."""

    epsilon: float
    sigma: float

    def __init__(self, element_types, epsilon=1.0, sigma=1.0, rcutfac=1.25):
        # ARGS: interface, element_types, ndescriptors, nparams, rcutfac
        super().__init__(None, element_types, 1, 3, rcutfac)
        # Mimicking the LJ pair-style:
        # pair_style lj/cut 2.5
        # pair_coeff * * 1 1
        self.epsilon = epsilon
        self.sigma = sigma

    def compute_gradients(self, data):
        """Test compute_gradients."""

    def compute_descriptors(self, data):
        """Test compute_descriptors."""

    def compute_forces(self, data):
        """Test compute_forces."""

        # NOTE: Use data.rij_max with JAX.
        #       dlpack requires cudnn:
        rij = jax.dlpack.from_dlpack(data.rij_max.toDlpack())
        eij, fij = lj_potential(self.epsilon, self.sigma, rij)

        # Convert back to cupy.
        eij = cupy.from_dlpack(jax.dlpack.to_dlpack(eij)).astype(np.float64)
        fij = cupy.from_dlpack(jax.dlpack.to_dlpack(fij)).astype(np.float64)

        # Send to LAMMPS.
        data.update_pair_energy(eij)
        data.update_pair_forces(fij)
