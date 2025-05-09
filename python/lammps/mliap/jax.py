from functools import partial

import jax
import jax.extend
import jax.numpy as jnp
from jax.interpreters import ad, batching, mlir, xla

forward_exchange_p = jax.extend.core.Primitive("forward_exchange")
reverse_exchange_p = jax.extend.core.Primitive("reverse_exchange")


def forward_exchange(x: jax.Array, pair_handle) -> jax.Array:
    return forward_exchange_p.bind(x, pair_handle=pair_handle)


def reverse_exchange(x: jax.Array, pair_handle) -> jax.Array:
    return reverse_exchange_p.bind(x, pair_handle=pair_handle)


forward_exchange_p.def_abstract_eval(lambda x, **_: x)
reverse_exchange_p.def_abstract_eval(lambda x, **_: x)
forward_exchange_p.def_impl(partial(xla.apply_primitive, forward_exchange_p))
reverse_exchange_p.def_impl(partial(xla.apply_primitive, reverse_exchange_p))


def forward_exchange_impl(x, pair_handle):
    if x.dtype == jax.numpy.float32:
        return jax.ffi.ffi_call("forward_exchange_fp32", x)(x, pair_handle=pair_handle)
    elif x.dtype == jax.numpy.float64:
        return jax.ffi.ffi_call("forward_exchange_fp64", x)(x, pair_handle=pair_handle)
    else:
        raise ValueError(f"Unsupported dtype: {x.dtype}")


mlir.register_lowering(
    forward_exchange_p, mlir.lower_fun(forward_exchange_impl, False), "cuda"
)


def reverse_exchange_impl(x, pair_handle):
    if x.dtype == jax.numpy.float32:
        return jax.ffi.ffi_call("reverse_exchange_fp32", x)(x, pair_handle=pair_handle)
    elif x.dtype == jax.numpy.float64:
        return jax.ffi.ffi_call("reverse_exchange_fp64", x)(x, pair_handle=pair_handle)
    else:
        raise ValueError(f"Unsupported dtype: {x.dtype}")


mlir.register_lowering(
    reverse_exchange_p, mlir.lower_fun(reverse_exchange_impl, False), "cuda"
)


# Automatic differentiation rules (JVP and transpose)
def jvp(
    primals: tuple[jax.Array, ...],
    tangents: tuple[jax.Array | ad.Zero, ...],
    pair_handle: int,
) -> tuple[jax.Array, jax.Array | ad.Zero]:
    (primal,) = primals
    (tangent,) = tangents
    out = forward_exchange_p.bind(primal, pair_handle=pair_handle)
    dout = forward_exchange_p.bind(tangent, pair_handle=pair_handle)
    return out, dout


ad.primitive_jvps[forward_exchange_p] = jvp


def jvp(
    primals: tuple[jax.Array, ...],
    tangents: tuple[jax.Array | ad.Zero, ...],
    pair_handle: int,
) -> tuple[jax.Array, jax.Array | ad.Zero]:
    (primal,) = primals
    (tangent,) = tangents
    out = reverse_exchange_p.bind(primal, pair_handle=pair_handle)
    dout = reverse_exchange_p.bind(tangent, pair_handle=pair_handle)
    return out, dout


ad.primitive_jvps[reverse_exchange_p] = jvp


def transpose(
    cotangent: jax.Array | ad.Zero,
    input: jax.Array | ad.UndefinedPrimal,
    pair_handle: int,
) -> tuple[jax.Array | ad.Zero, ...]:
    assert ad.is_undefined_primal(input)
    return (reverse_exchange_p.bind(cotangent, pair_handle=pair_handle),)


ad.primitive_transposes[forward_exchange_p] = transpose


def transpose(
    cotangent: jax.Array | ad.Zero,
    input: jax.Array | ad.UndefinedPrimal,
    pair_handle: int,
) -> tuple[jax.Array | ad.Zero, ...]:
    assert ad.is_undefined_primal(input)
    return (forward_exchange_p.bind(cotangent, pair_handle=pair_handle),)


ad.primitive_transposes[reverse_exchange_p] = transpose


def batch(
    primitive: jax.extend.core.Primitive,
    inputs: tuple[jax.Array, ...],
    batch_axes: tuple[int | None, ...],
    pair_handle: int,
) -> tuple[tuple[jax.Array, ...], tuple[int, ...]]:
    assert isinstance(inputs, tuple)
    assert isinstance(batch_axes, tuple)
    assert isinstance(pair_handle, int)
    (input,) = inputs
    (axis,) = batch_axes
    assert axis is not None
    input = jnp.moveaxis(input, axis, 1)
    out = primitive.bind(input, pair_handle=pair_handle)
    return out, 1


batching.primitive_batchers[forward_exchange_p] = partial(batch, forward_exchange_p)
batching.primitive_batchers[reverse_exchange_p] = partial(batch, reverse_exchange_p)
