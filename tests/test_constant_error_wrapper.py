import pytest


def test_constant_not_found_error_has_helpful_message(tmp_path):
    """Requesting a missing constant after loading a sparse profile raises
    ConstantNotFoundError with profile name and remediation hints."""
    from climt import (
        ConstantNotFoundError,
        get_constant_checked,
        load_atmospheric_properties,
        reset_atmospheric_properties,
    )

    toml = """
[planetary]
gravitational_acceleration = { value = 3.721, units = "m/s^2" }
"""
    p = tmp_path / "sparse.toml"
    p.write_text(toml)
    load_atmospheric_properties(str(p))

    try:
        with pytest.raises(ConstantNotFoundError) as exc_info:
            get_constant_checked("molar_mass_of_water_vapor", "g/mol")
        msg = str(exc_info.value)
        assert "molar_mass_of_water_vapor" in msg
        assert "sparse" in msg or str(p) in msg
        assert "set_constant" in msg or "profile" in msg.lower()
    finally:
        reset_atmospheric_properties()
