import numpy as np
import textwrap


def _write_fake_cia(tmp_path, name="N2-N2_fake.cia"):
    """Write a 2-block, 3-point synthetic HITRAN CIA file."""
    text = textwrap.dedent("""\
        N2-N2          10.0     30.0    3  200.0  1.000E-50    Karman2019
            10.0      1.0E-46
            20.0      2.0E-46
            30.0      4.0E-46
        N2-N2          10.0     30.0    3  300.0  1.000E-50    Karman2019
            10.0      2.0E-46
            20.0      4.0E-46
            30.0      8.0E-46
        """)
    path = tmp_path / name
    path.write_text(text)
    return path


def test_load_cia_blocks(tmp_path):
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.cork_table_builder.cia import load_cia_blocks

    p = _write_fake_cia(tmp_path)
    blocks = load_cia_blocks(str(p))
    assert sorted(blocks.keys()) == [200.0, 300.0]
    nu, k = blocks[200.0]
    np.testing.assert_allclose(nu, [10.0, 20.0, 30.0])
    np.testing.assert_allclose(k, [1e-46, 2e-46, 4e-46])


def test_cia_kappa_on_grid(tmp_path):
    """cia_kappa_on_grid returns kappa(T,p,nu) in m^2/kg for the pair."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.cork_table_builder.cia import cia_kappa_on_grid

    p = _write_fake_cia(tmp_path)
    T_grid = np.array([200.0, 250.0, 300.0])
    p_grid = np.array([1e4, 1e5])
    nu_grid = np.array([15.0, 25.0])
    kappa = cia_kappa_on_grid(
        str(p), pair=("N2", "N2"),
        vmr_a=1.0, vmr_b=1.0,
        background_gas="N2",
        absorbers={"N2": 1.0},  # used for mass-fraction normalisation
        T_grid=T_grid, p_grid=p_grid, nu_grid=nu_grid,
    )
    assert kappa.shape == (3, 2, 2)
    # CIA grows like p² at fixed T (amagat² scaling), so high-p > low-p
    assert (kappa[:, 1, :] > kappa[:, 0, :]).all()
    # And grows with T in the synthetic file
    assert kappa[2, 0, 0] > kappa[0, 0, 0]
