import numpy as np
import pandas as pd
import xarray as xr
from glob import glob

import sys
from config import *

import filter_mode

files = sorted(glob(str(OLR_DIR) + "/*.nc"))
var_name = "txuptp"
spd = 8

fileo = "/scratch/gb02/mr4682/data/txuptp/clim/txuptp.clim.nc"

def calc_clim_olr(files, var_name, spd):
    if len(files) == 1:
        ds = xr.open_dataset(files[0])
    elif len(files) > 1:
        ds = xr.open_mfdataset(files)
    else:
        sys.exit("No .nc files are found!")
    
    olr = ds[var_name]
    
    print(olr)
    
    if len(olr.dims) == 3:
        clim = filter_mode.clim.calcClimTLL(olr, spd=spd)
        clim = filter_mode.clim.smthClimTLL(clim, spd=spd, nsmth=4)
    elif len(olr.dims) == 4:
        clim = filter_mode.clim.calcClimTLLL(olr, spd=spd)
        clim = filter_mode.clim.smthClimTLLL(clim, spd=spd, nsmth=4)
    else:
        sys.exit("Only 3D or 4D arrays are accepted!")
    
    print(clim)
    
    return clim

if __name__ == "__main__":
    clim = calc_clim_olr(files, var_name, spd)
    clim.to_netcdf(path=fileo, format="NETCDF4")
