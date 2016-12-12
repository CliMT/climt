from .units import to_units as to_units_function
import xarray as xr


class DataArray(xr.DataArray):

    def __add__(self, other):
        result = super(DataArray, self).__add__(other)
        result.attrs = self.attrs
        return result

    def __sub__(self, other):
        result = super(DataArray, self).__sub__(other)
        result.attrs = self.attrs
        return result

    def to_units(self, units):
        return to_units_function(self, units)
