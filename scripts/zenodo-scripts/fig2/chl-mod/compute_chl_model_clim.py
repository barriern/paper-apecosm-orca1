import os.path
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from glob import glob
import os

dirin = os.path.join('..', '..', 'data', 'nemo-pisces')
dirin

data = xr.open_mfdataset(os.path.join(dirin, '*add_T*nc'), combine="by_coords")
data = data.isel(olevel=0)
data

chl = data['CHL']
chl = chl.sel(time_counter=slice('1998-01-01', '2018-12-31'))
chl

clim = chl.groupby('time_counter.month').mean(dim='time_counter')
clim.name = 'clim_chl'

# +
delayed = clim.to_netcdf('clim_chl_monthly_model.nc', compute=False)

from dask.diagnostics import ProgressBar
with ProgressBar():
    delayed.compute()