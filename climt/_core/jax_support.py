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


def require_x64() -> None:
    """Raise JaxPrecisionError unless jax.config.jax_enable_x64 is True. The
    message contains: jax.config.update("jax_enable_x64", True)"""
    import jax

    if not jax.config.jax_enable_x64:
        raise JaxPrecisionError(
            "climt's JAX kernels require 64-bit floating-point precision; enable it "
            'before creating JAX arrays with jax.config.update("jax_enable_x64", True)'
        )
