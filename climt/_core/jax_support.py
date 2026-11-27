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


def count_jaxpr_equations(function, *args) -> int:
    """Trace function(*args) with jax.make_jaxpr and return the number of
    equations, descending into every nested jaxpr (jit, scan, cond, while
    bodies) and counting each nested jaxpr object once."""
    import jax
    try:
        from jax.extend.core import ClosedJaxpr, Jaxpr
    except ImportError:  # pragma: no cover - older JAX fallback
        from jax.core import ClosedJaxpr, Jaxpr

    closed = jax.make_jaxpr(function)(*args)
    seen = set()

    def _iter_nested(obj):
        if isinstance(obj, ClosedJaxpr):
            yield obj.jaxpr
        elif isinstance(obj, Jaxpr):
            yield obj
        elif isinstance(obj, (tuple, list)):
            for item in obj:
                yield from _iter_nested(item)
        elif isinstance(obj, dict):
            for val in obj.values():
                yield from _iter_nested(val)

    def _count(jaxpr) -> int:
        key = id(jaxpr)
        if key in seen:
            return 0
        seen.add(key)
        total = len(jaxpr.eqns)
        for eqn in jaxpr.eqns:
            for val in eqn.params.values():
                for sub in _iter_nested(val):
                    total += _count(sub)
        return total

    return _count(closed.jaxpr)
