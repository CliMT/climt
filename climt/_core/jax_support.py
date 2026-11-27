"""Helpers for optional JAX dispatch and precision checks.

This module must not import ``jax`` at module level so that ``import climt``
stays JAX-free whether or not JAX is installed.
"""

import sys


class JaxPrecisionError(RuntimeError):
    """Raised when a JAX kernel is called with 64-bit precision disabled."""


def is_jax_array(value) -> bool:
    """True for a jax.Array (including tracers). Never imports jax: if 'jax'
    is not in sys.modules the answer is False."""
    jax = sys.modules.get("jax")
    if jax is None:
        return False
    return isinstance(value, jax.Array)


def require_x64() -> None:
    """Raise JaxPrecisionError unless jax.config.jax_enable_x64 is True. The
    message contains: jax.config.update("jax_enable_x64", True)"""
    import jax

    if not jax.config.jax_enable_x64:
        raise JaxPrecisionError(
            "CliMT's JAX kernels require 64-bit floating-point precision, "
            "but JAX x64 mode is disabled. Enable it before creating JAX "
            'arrays with: jax.config.update("jax_enable_x64", True)'
        )
