import nctools as nc
import sys

lon=0.0
lat=0.0
if __name__ == '__main__':
    if (len(sys.argv) > 1):
        szen = float(sys.argv[1])
    else:
        raise RuntimeError('please enter a zenith angle')

nc.ncout2d('case6.szen', lon, lat, szen, longname = 'Solar Zenith Angle', units = 'degrees')
