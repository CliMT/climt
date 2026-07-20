import numpy as np
from climt import get_grid


def test_soil_grid_present():
    grid = get_grid(nx=4, ny=2, nz=10, n_soil_interface_levels=4)
    assert "height_on_soil_interface_levels" in grid
    assert grid["height_on_soil_interface_levels"].shape[0] == 4
