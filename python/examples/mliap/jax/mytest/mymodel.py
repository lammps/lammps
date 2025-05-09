import jax
import jax.numpy as jnp
from lammps.mliap.jax import forward_exchange, reverse_exchange
from lammps.mliap.mliap_unified_abc import MLIAPUnified
from functools import partial


class MLIAPMod(MLIAPUnified):
    def __init__(self, element_types=None):
        super().__init__()
        self.ndescriptors = 1
        self.element_types = element_types
        self.rcutfac = 1.0  # Half of radial cutoff
        self.nparams = 1

    def compute_gradients(self, data):
        pass

    def compute_descriptors(self, data):
        pass

    def compute_forces(self, data):
        natoms = data.nlocal
        nghosts = data.ntotal - data.nlocal

        # call once this
        data.register_jax_ffi_exchange()
        pair_handle = data.get_pair_handle()

        @partial(jax.jit, static_argnums=(1,))
        def f(node_features, pair_handle):
            return forward_exchange(node_features, pair_handle)

        node_features = jnp.concatenate(
            [jnp.ones((natoms, 4)), jnp.zeros((nghosts, 4))],
            axis=0,
        )
        node_features = f(node_features, pair_handle)
        print("forward", node_features[0], node_features[-1])

        grad = jax.vjp(lambda x: f(x, pair_handle), node_features)[1](node_features)
        print("grad", grad[0][0], grad[0][-1])

        node_features = jnp.concatenate(
            [jnp.ones((2, natoms, 4)), jnp.zeros((2, nghosts, 4))],
            axis=1,
        )
        node_features = jax.vmap(lambda x: f(x, pair_handle))(node_features)
        print("vmap", node_features[:, 0], node_features[:, -1])
