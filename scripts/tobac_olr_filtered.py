import numpy as np
import pandas as pd
import xarray as xr
import iris
import tobac
import os
from glob import glob

import sys

import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)
warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)
warnings.filterwarnings('ignore', category=FutureWarning, append=True)
warnings.filterwarnings('ignore',category=pd.io.pytables.PerformanceWarning)

mode = "IG Wave"

diri = "/scratch/gb02/mr4682/data/regridded/UM/not_removing_tc/"

if mode == "Moisture Mode":
    mode_str = "moisture.mode"
    mode_var_name = "moisture_mode"
    v_max = 9.0
elif mode == "Mixed System":
    mode_str = "mixed.system"
    mode_var_name = "mixed_system"
    v_max = 27.0
elif mode == "IG Wave":
    mode_str = "ig.wave"
    mode_var_name = "ig_wave"
    v_max = 50.0
elif mode == "Tropical Cyclone":
    modes_str = "tc"
    mode_var_name = mode_str
    v_max = 10.0
else:
    sys.exit(f"There is no mode named {mode}!")

filei = f"olr.zoom10.to.0p25deg.{mode_str}"

diro = "/scratch/gb02/mr4682/data/regridded/UM/tracks/"
fileo = f"Track.{filei}.csv"

time_start = "2020-01-20 00:00:00"
time_end = "2021-04-01 00:00:00"

ds = xr.open_dataset(f"{diri}{filei}.nc")

time = ds["time"].sel(time=slice(time_start, time_end))
latitude = ds["latitude"]
longitude = ds["longitude"]

olr = ds[f"rlut_{mode_var_name}"].sel(time=slice(time_start, time_end)).compute()

olr_stddev = olr.std(ddof=1)

olr_standardised = olr / olr_stddev
olr_standardised = olr_standardised.fillna(0.0)

print(olr_standardised)

dxy = np.radians(longitude.values[1] - longitude.values[0]) * 6.371e+06
dt = (time.values[1] - time.values[0]) / np.timedelta64(1, "s")

parameters_features = {}
parameters_features["threshold"] = [-2.0, -2.5, -3.0, -3.5, -4.0]
parameters_features["target"] = "minimum"
parameters_features["position_threshold"] = "weighted_diff"
parameters_features["sigma_threshold"] = 1.0
parameters_features["n_min_threshold"] = 4
parameters_features["PBC_flag"] = "hdim_2"

Features = tobac.feature_detection_multithreshold(olr_standardised, dxy, **parameters_features)

parameters_linking = {}
parameters_linking["v_max"] = v_max
parameters_linking["stubs"] = 3
parameters_linking["method_linking"] = "predict"
parameters_linking["PBC_flag"] = "hdim_2"
parameters_linking["min_h2"] = 0
parameters_linking["max_h2"] = len(longitude) - 1

Track = tobac.linking_trackpy(Features, olr_standardised, dt=dt, dxy=dxy, **parameters_linking)
Track = Track.loc[Track["cell"] != -1]

print(Track)

Track.to_csv(f"{diro}{fileo}")
