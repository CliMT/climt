"""Surface albedo process (BEST Eqs 5.5-5.8)."""


class SurfaceAlbedo:
    def __call__(self, soil_props, wetness, area_type):
        raise NotImplementedError


class BestSurfaceAlbedo(SurfaceAlbedo):
    def __call__(self, soil_props, wetness, area_type):
        if area_type == "land_ice":
            alpha_sw = 0.60 + 0.06 * (1.0 - wetness)
            alpha_lw = alpha_sw / 3.0
        else:
            alpha_sw = 0.10 + 0.1 * soil_props["colour"] + 0.06 * (1.0 - wetness)
            alpha_lw = 2.0 * alpha_sw
        return {"alpha_sw": alpha_sw, "alpha_lw": alpha_lw}
