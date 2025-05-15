import numpy as np
import pandas as pd
import xarray as xr
import filter_mode

import sys

mode = "IG Wave"
margin = 2.0

diri = "/scratch/gb02/mr4682/data/regridded/UM/not_removing_tc/"
diro = "/scratch/gb02/mr4682/data/regridded/UM/tracks/"

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

filei = f"olr.zoom10.to.0p25deg.{mode_str}"

ds = xr.open_dataset(f"{diri}{filei}.nc")

time = ds["time"]
latitude = ds["latitude"]
longitude = ds["longitude"]

olr = ds[f"rlut_{mode_var_name}"].compute()

print(olr)

Track = pd.read_csv(f"{diro}Track.{filei}.csv")

print(Track)

max_activity = filter_mode.track.calc_max_activity(olr, Track, latitude, longitude, margin)

print(max_activity)

max_activity.to_csv(f"{diro}max.activity.{filei}.csv")
