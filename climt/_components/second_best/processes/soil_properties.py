"""Soil property process (BEST Eqs 4.10-4.12)."""


class SoilProperties:
    """Contract: __call__(soil_type, area_type) -> dict of soil parameters."""
    def __call__(self, soil_type, area_type):
        raise NotImplementedError


class BestSoilProperties(SoilProperties):
    _COLOUR = {"clay": 0.2, "sand": 1.0}
    _TEXTURE = {"clay": 0.0, "sand": 9.0}
    _B = {"clay": 10.0, "sand": 4.0}
    _K_H0 = {"clay": 0.001, "sand": 0.1}

    def __call__(self, soil_type, area_type):
        colour = self._COLOUR[soil_type]
        texture = self._TEXTURE[soil_type]
        if area_type == "land_ice":
            texture = 0.07
        porosity = 0.6 - 0.03 * texture
        field_capacity = (0.95 - 0.086 * texture) * porosity
        wilting_point = 0.01 if area_type == "land_ice" else porosity - 0.03
        return {
            "colour": colour, "texture": texture, "porosity": porosity,
            "field_capacity": field_capacity, "wilting_point": wilting_point,
            "B": self._B[soil_type], "K_H0": self._K_H0[soil_type],
            "psi_0": -0.2,
        }
