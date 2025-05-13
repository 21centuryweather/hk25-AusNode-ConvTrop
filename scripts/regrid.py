import numpy as np
import xarray as xr
import healpy as hp

res = 0.25
res_str = "0p25"

zoom = '10'
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

olr = ds2d["rlut"]

nside = hp.get_nside(olr)

lon = np.arange(-180.0, 180.0, res)
lat = np.arange(-30.0, 30.0 + res, res)

cells = get_nn_lon_lat_index(nside, lon, lat)

olr_regridded = olr.isel(cell=cells)
olr_regridded = olr_regridded.rename({"lon": "longitude", "lat": "latitude"})

del olr_regridded.attrs["hiopy::time_method"]
del olr_regridded.attrs["hiopy::nnn"]
del olr_regridded.attrs["hiopy::enable"]

fileo = f"/scratch/gb02/mr4682/data/regridded/ICON/olr.zoom{zoom}.to.{res_str}deg.nc"

olr_regridded.to_netcdf(path=fileo)