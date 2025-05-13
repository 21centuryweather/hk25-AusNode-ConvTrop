import numpy as np
import pandas as pd
import xarray as xr
from glob import glob

import sys
from config import *

import filter_mode

files = sorted(glob("/scratch/gb02/mr4682/data/regridded/ICON/olr.anom.zoom10.to.0p25deg.nc"))
var_name = "rlut"
spd = 8

diro = "/scratch/gb02/mr4682/data/regridded/ICON/not_removing_tc/"

def filter_mode_olr(files, var_name, spd, mode):
    if len(files) == 1:
        ds = xr.open_dataset(files[0])
    elif len(files) > 1:
        ds = xr.open_mfdataset(files)
    else:
        sys.exit("No .nc files are found!")

    olr = ds[f"{var_name}_anom"]

    if len(olr.dims) == 3:
        olr = filter_mode.filter.pad_zeros_TLL(olr, spd)
        olr_filtered = filter_mode.filter.filter_mode_TLL(olr, mode, len(olr["time"].values), spd)
    elif len(olr.dims) == 4:
        olr = filter_mode.filter.pad_zeros_TLLL(olr, spd)
        olr_filtered = filter_mode.filter.filter_mode_TLLL(olr, mode, len(olr["time"].values), spd)
    else:
        sys.exit("Only 3D or 4D arrays are accepted!")

    print(olr_filtered)

    return olr_filtered

if __name__ == "__main__":
    olr_moisture_mode = filter_mode_olr(files, var_name, spd, "Moisture Mode")
    olr_mixed_system = filter_mode_olr(files, var_name, spd, "Mixed System")
    olr_ig_wave = filter_mode_olr(files, var_name, spd, "IG Wave")

    dso1 = xr.Dataset()
    dso1[f"{var_name}_moisture_mode"] = olr_moisture_mode
    
    dso2 = xr.Dataset()
    dso2[f"{var_name}_mixed_system"] = olr_mixed_system
    
    dso3 = xr.Dataset()
    dso3[f"{var_name}_ig_wave"] = olr_ig_wave

    dso1.to_netcdf(path=f"{diro}olr.moisture.mode.zoom10.to.0p25deg.nc", format="NETCDF4")
    dso2.to_netcdf(path=f"{diro}olr.mixed.system.zoom10.to.0p25deg.nc", format="NETCDF4")
    dso3.to_netcdf(path=f"{diro}olr.ig.wave.zoom10.to.0p25deg.nc", format="NETCDF4")
