import numpy as np
import pandas as pd
import xarray as xr
import filter_mode

import sys

mode = "Moisture Mode"

diri = "/scratch/gb02/mr4682/data/regridded/UM/"
diro = "/scratch/gb02/mr4682/data/regridded/UM/composites/"

filei = f"{diri}olr.zoom10.to.0p25deg.anom.nc"
var_name = "rlut_anom"

ds = xr.open_dataset(filei)

time_start = "2020-01-20 00:00:00"
time_end = "2021-04-01 00:00:00"

time = ds["time"].sel(time=slice(time_start, time_end))
latitude = ds["latitude"]
longitude = ds["longitude"]

field = ds[var_name].sel(time=slice(time_start, time_end))

if mode == "Moisture Mode":
    mode_str = "moisture.mode"
    mode_var_name = "moisture_mode"
elif mode == "Mixed System":
    mode_str = "mixed.system"
    mode_var_name = "mixed_system"
elif mode == "IG Wave":
    mode_str = "ig.wave"
    mode_var_name = "ig_wave"
elif mode == "Tropical Cyclone":
    modes_str = "tc"
    mode_var_name = mode_str
else:
    sys.exit(f"There is no mode named {mode}!")

track = pd.read_csv(f"{diri}tracks/Track.olr.zoom10.to.0p25deg.{mode_str}.csv")

print(track)

max_activity = pd.read_csv(f"{diri}tracks/max.activity.olr.zoom10.to.0p25deg.{mode_str}.csv")

print(max_activity)

cell_to_use = []

for cell in max_activity["cell"].values:
    track_i = track.loc[track["cell"] == cell]
    if ~(np.any(np.isnan(track_i["latitude"].values)) | np.any(np.isnan(track_i["longitude"].values))):
        cell_to_use.append(cell)

track = track.loc[np.isin(track["cell"], cell_to_use)]

print(track)

max_activity = max_activity.loc[np.isin(max_activity["cell"], cell_to_use)]

print(max_activity)

season = "ALL"  # ["ALL", "DJF", "MAM", "JJA", "SON"]
domain = "EQ"   # ["EQ", "NH", "SH"]

x_size = 15.0
y_size = 15.0

grid_res = longitude.values[1] - longitude.values[0]
time_res = (time.values[1] - time.values[0]) / np.timedelta64(1, "h")

lagmax = 40

if len(field.dims) == 3:
    composite, pvalue = filter_mode.track.compositeTLL(field, track, max_activity, latitude, longitude, season, domain, x_size, y_size, grid_res, lagmax, time_res)
elif len(field.dims) == 4:
    composite, pvalue = filter_mode.track.compositeTLLL(field, track, max_activity, latitude, longitude, season, domain, x_size, y_size, grid_res, lagmax, time_res)

print(composite)

dso = xr.Dataset()
dso[f"composite_{var_name}"] = composite
dso[f"pvalue_{var_name}"] = pvalue

fileo = f"{diro}composite_{var_name}_{mode_var_name}_{season}_{domain}.nc"

dso.to_netcdf(fileo, format="NETCDF4")
