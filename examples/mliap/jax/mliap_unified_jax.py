from lammps.mliap.mliap_unified_abc import MLIAPUnified
import numpy as np
import jax
import jax.numpy as jnp
from jax import jit
from functools import partial
import os

# Required else get `jaxlib.xla_extension.XlaRuntimeError: RESOURCE_EXHAUSTED: Out of memory`
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"]="false"
os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"]=".XX"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"]="platform"

@jax.jit
def lj_potential(epsilon, sigma, rij):
    def _tot_e(rij):
        """A differentiable fn for total energy."""
        r2inv = 1.0 / jnp.sum(rij ** 2, axis=1)
        r6inv = r2inv * r2inv * r2inv

        lj1 = 4.0 * epsilon * sigma**12
        lj2 = 4.0 * epsilon * sigma**6

        eij = r6inv * (lj1 * r6inv - lj2)
        return 0.5 * jnp.sum(eij), eij
    # Compute _tot_e and its derivative.
    (_, eij), fij = jax.value_and_grad(_tot_e, has_aux=True)(rij)
    return eij, fij


class MLIAPUnifiedJAX(MLIAPUnified):
    """Test implementation for MLIAPUnified."""

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
        rij = data.rij_max

        eij, fij = lj_potential(self.epsilon, self.sigma, rij)

        data.update_pair_energy(np.array(eij, dtype=np.float64))
        data.update_pair_forces(np.array(fij, dtype=np.float64))
