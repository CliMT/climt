import pytest

jax = pytest.importorskip("jax")
pytestmark = pytest.mark.jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")

import jax.numpy as jnp                                  # noqa: E402
from climt._core.jax_support import count_jaxpr_equations   # noqa: E402


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
    inner = jax.jit(lambda x: jnp.sin(x) + jnp.cos(x))     # sin, cos, add -> 3
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
