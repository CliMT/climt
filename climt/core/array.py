try:
    import xarray as xr
except ImportError:
    xr = None

if xr is not None:
    class DataArray(xr.DataArray):

        def __add__(self, other):
            result = super(DataArray, self).__add__(other)
            result.attrs = self.attrs
            return result

        def __sub__(self, other):
            result = super(DataArray, self).__sub__(other)
            result.attrs = self.attrs
            return result
