import numpy as np
import xarray as xr
import healpy as hp

res = 0.25
res_str = "0p25"

variable = "ua"
pressure = 850

zoom = '9'
file = '/g/data/qx55/germany_node/d3hp003.zarr/PT3H_mean_z' + zoom + '_atm.zarr'

ds2d = xr.open_zarr(file)

def get_nn_lon_lat_index(nside, lons, lats):
    """
    nside: integer, power of 2. The return of hp.get_nside()
    lons: uniques values of longitudes
    lats: uniques values of latitudes
    returns: array with the HEALPix cells that are closest to the lon/lat grid
    """
    lons2, lats2 = np.meshgrid(lons, lats)
    return xr.DataArray(
        hp.ang2pix(nside, lons2, lats2, nest = True, lonlat = True),
        coords=[("lat", lats), ("lon", lons)],
    )

x = ds2d[variable].sel(pressure=pressure)

nside = hp.get_nside(x)

lon = np.arange(-180.0, 180.0, res)
lat = np.arange(-30.0, 30.0 + res, res)

cells = get_nn_lon_lat_index(nside, lon, lat)

x_regridded = x.isel(cell=cells)
x_regridded = x_regridded.rename({"lon": "longitude", "lat": "latitude"})

fileo = f"/scratch/gb02/mr4682/data/regridded/UM/{variable}.{pressure}.zoom{zoom}.to.{res_str}deg.nc"

x_regridded.to_netcdf(path=fileo)
