import pytest
import numpy as np


def test_load_earth_profile():
    """load_atmospheric_properties('earth') sets Earth constants."""
    from climt import load_atmospheric_properties
    from sympl import get_constant

    load_atmospheric_properties("earth")
    g = get_constant("gravitational_acceleration", "m/s^2")
    assert abs(g - 9.80665) < 1e-3


def test_load_mars_profile():
    """load_atmospheric_properties('mars') sets Mars constants."""
    from climt import load_atmospheric_properties
    from sympl import get_constant

    load_atmospheric_properties("mars")
    g = get_constant("gravitational_acceleration", "m/s^2")
    assert abs(g - 3.721) < 0.01
    cp = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/K")
    assert abs(cp - 735.0) < 10.0


def test_reset_restores_previous():
    """reset_atmospheric_properties restores prior constants."""
    from climt import load_atmospheric_properties, reset_atmospheric_properties
    from sympl import get_constant

    load_atmospheric_properties("earth")
    g_earth = get_constant("gravitational_acceleration", "m/s^2")

    load_atmospheric_properties("mars")
    g_mars = get_constant("gravitational_acceleration", "m/s^2")
    assert abs(g_mars - 3.721) < 0.01

    reset_atmospheric_properties()
    g_restored = get_constant("gravitational_acceleration", "m/s^2")
    assert abs(g_restored - g_earth) < 1e-6


def test_load_custom_toml(tmp_path):
    """load_atmospheric_properties with a file path loads custom TOML."""
    from climt import load_atmospheric_properties
    from sympl import get_constant

    toml_content = """
[planetary]
gravitational_acceleration = { value = 24.79, units = "m/s^2" }

[bulk_atmosphere]
molar_mass_of_dry_air = { value = 2.2, units = "g/mol" }
"""
    p = tmp_path / "jupiter.toml"
    p.write_text(toml_content)

    load_atmospheric_properties(str(p))
    g = get_constant("gravitational_acceleration", "m/s^2")
    assert abs(g - 24.79) < 0.01


def test_load_nonexistent_profile_raises():
    """Unknown profile name raises FileNotFoundError."""
    from climt import load_atmospheric_properties

    with pytest.raises(FileNotFoundError):
        load_atmospheric_properties("nonexistent_planet_xyz")
