"""Compile-time and traced-graph-size budget tests for the Emanuel JAX kernel.

Measured on ``arm64`` (Apple Silicon CPU) with JAX ``0.10.1`` at ``nlev=30, ncol=1``:
- Baseline (as merged, Task 3):
  - Forward:  63,262 equations, 16.392 s first call (0.624 ms steady median)
  - Gradient: 180,701 equations, 719.943 s first call (4.834 ms steady median)
- Achieved (after loop regrouping, Task 4):
  - Forward:  3,123 equations, 0.751 s first call (0.108 ms steady median)
  - Gradient: 5,392 equations, 2.122 s first call (0.337 ms steady median)

Pinned limits are 1.10x the achieved equation counts rounded up to the next 100
(3,500 forward; 6,000 gradient). Neither graph-size test exceeds 1 s, so both
run in the default (non-slow) ``jax`` suite.
"""
import os
import time

import pytest

jax = pytest.importorskip("jax")
pytestmark = pytest.mark.jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")

import jax.numpy as jnp                                  # noqa: E402
from climt._core.jax_support import count_jaxpr_equations   # noqa: E402
from tests.jax_cases import emanuel_case, emanuel_forward, emanuel_loss  # noqa: E402

FORWARD_EQUATION_LIMIT = 3500
GRADIENT_EQUATION_LIMIT = 6000
SCALE = float(os.environ.get("CLIMT_JAX_BUDGET_SCALE", "1.0"))


def test_count_is_exact_on_a_known_program():
    def body(carry, x):
        return carry + x * 2.0, carry

    def program(xs):
        return jax.lax.scan(body, 0.0, xs)[0] + 1.0

    # jax.make_jaxpr(program)(jnp.arange(5.0)) on JAX 0.10.1:
    # { lambda ; a:f64[5]. let
    #     b:f64[] _:f64[5] = scan[
    #       jaxpr={ lambda ; c:f64[] d:f64[]. let
    #           e:f64[] = mul d 2.0:f64[]
    #           f:f64[] = convert_element_type[new_dtype=float64 weak_type=False] c
    #           g:f64[] = add f e
    #         in (g, c) }
    #       length=5 num_carry=1 num_consts=0 reverse=False unroll=1
    #     ] 0.0:f64[] a
    #     h:f64[] = add b 1.0:f64[]
    #   in (h,) }
    # top level: scan, add -> 2; scan body: mul, convert_element_type, add -> 3
    assert count_jaxpr_equations(program, jnp.arange(5.0)) == 5


def test_count_descends_into_jitted_functions():
    inner = jax.jit(lambda x: jnp.sin(x) + jnp.cos(x))
    # jax.make_jaxpr(lambda x: inner(x) * 2.0)(1.0) on JAX 0.10.1:
    # { lambda ; a:f64[]. let
    #     b:f64[] = jit[
    #       name=<lambda>
    #       jaxpr={ lambda ; a:f64[]. let
    #           c:f64[] = sin a
    #           d:f64[] = cos a
    #           b:f64[] = add c d
    #         in (b,) }
    #     ] a
    #     e:f64[] = mul b 2.0:f64[]
    #   in (e,) }
    # top level: jit, mul -> 2; inner jit body: sin, cos, add -> 3
    assert count_jaxpr_equations(lambda x: inner(x) * 2.0, 1.0) == 5


def test_count_descends_into_cond_while_and_deduplicates_nested_jaxpr():
    inner = jax.jit(lambda x: jnp.sin(x) + jnp.cos(x))  # 3 eqns inside

    def program(x):
        y = inner(x) + inner(x)
        y = jax.lax.cond(y > 0.0, lambda z: z * 2.0, lambda z: z - 1.0, y)
        return jax.lax.while_loop(lambda z: z < 10.0, lambda z: z + 1.0, y)

    closed = jax.make_jaxpr(program)(jnp.asarray(1.0))
    seen = set()

    def manual(j):
        if id(j) in seen:
            return 0
        seen.add(id(j))
        total = len(j.eqns)
        for e in j.eqns:
            for v in e.params.values():
                items = v if isinstance(v, (tuple, list)) else (v,)
                for item in items:
                    sub = getattr(item, "jaxpr", item)
                    if hasattr(sub, "eqns"):
                        total += manual(sub)
        return total

    assert count_jaxpr_equations(program, jnp.asarray(1.0)) == manual(closed.jaxpr)


def test_forward_graph_size():
    component, state, timestep = emanuel_case()
    count = count_jaxpr_equations(emanuel_forward(component, state, timestep), state["air_temperature"])
    assert count <= FORWARD_EQUATION_LIMIT, count


def test_gradient_graph_size():
    component, state, timestep = emanuel_case()
    count = count_jaxpr_equations(jax.grad(emanuel_loss(component, state, timestep)), state["air_temperature"])
    assert count <= GRADIENT_EQUATION_LIMIT, count


@pytest.mark.slow
def test_forward_first_call_within_budget():
    component, state, timestep = emanuel_case()
    start = time.perf_counter()
    jax.block_until_ready(jax.jit(emanuel_forward(component, state, timestep))(state["air_temperature"]))
    assert time.perf_counter() - start <= 30.0 * SCALE


@pytest.mark.slow
def test_gradient_first_call_within_budget():
    component, state, timestep = emanuel_case()
    start = time.perf_counter()
    jax.block_until_ready(jax.jit(jax.grad(emanuel_loss(component, state, timestep)))(state["air_temperature"]))
    assert time.perf_counter() - start <= 120.0 * SCALE
