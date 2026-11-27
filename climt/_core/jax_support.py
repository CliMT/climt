# -*- coding: utf-8 -*-
import sys


class JaxPrecisionError(RuntimeError):
    """Raised when a climt JAX kernel is called without 64-bit precision enabled."""


def is_jax_array(value) -> bool:
    """True for a jax.Array (including tracers). Never imports jax: if 'jax'
    is not in sys.modules the answer is False."""
    jax_module = sys.modules.get("jax")
    if jax_module is None:
        return False
    array_type = getattr(jax_module, "Array", None)
    return array_type is not None and isinstance(value, array_type)


def require_x64(*arrays) -> None:
    """Raise JaxPrecisionError unless jax.config.jax_enable_x64 is True and
    any supplied floating-point JAX arrays are float64. The message contains:
    jax.config.update("jax_enable_x64", True)"""
    import jax
    import jax.numpy as jnp

    if not jax.config.jax_enable_x64:
        raise JaxPrecisionError(
            "climt's JAX kernels require 64-bit floating-point precision; enable it "
            'before creating JAX arrays with jax.config.update("jax_enable_x64", True)'
        )
    for arr in arrays:
        if arr is not None and is_jax_array(arr):
            dtype = getattr(arr, "dtype", None)
            if dtype is not None and jnp.issubdtype(dtype, jnp.floating) and dtype != jnp.float64:
                raise JaxPrecisionError(
                    f"climt's JAX kernels require float64 arrays, got {dtype}; enable "
                    "64-bit mode before creating JAX arrays with "
                    'jax.config.update("jax_enable_x64", True)'
                )


def count_jaxpr_equations(function, *args) -> int:
    """Trace function(*args) with jax.make_jaxpr and return the number of
    equations, descending into every nested jaxpr (jit, scan, cond, while
    bodies) and counting each nested jaxpr object once."""
    import jax

    closed = jax.make_jaxpr(function)(*args)
    seen_ids = set()
    seen_objs = []

    def _walk(obj) -> int:
        if hasattr(obj, "jaxpr") and hasattr(obj, "consts"):
            obj = obj.jaxpr
        if hasattr(obj, "eqns") and hasattr(obj, "invars"):
            if id(obj) in seen_ids:
                return 0
            seen_ids.add(id(obj))
            seen_objs.append(obj)
            total = len(obj.eqns)
            for eqn in obj.eqns:
                for param_val in eqn.params.values():
                    total += _walk(param_val)
            return total
        if isinstance(obj, (tuple, list)):
            return sum(_walk(item) for item in obj)
        if isinstance(obj, dict):
            return sum(_walk(item) for item in obj.values())
        return 0

    return _walk(closed)
