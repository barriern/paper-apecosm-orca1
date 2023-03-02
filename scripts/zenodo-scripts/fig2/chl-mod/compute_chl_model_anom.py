import os.path
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from glob import glob
import os

dirin = os.path.join('..', '..', 'data', 'nemo-pisces')

datac = xr.open_dataset('clim_chl_monthly_model.nc')
clim = datac['clim_chl']
clim

data = xr.open_mfdataset(os.path.join(dirin, '*add_T*nc'), combine="by_coords")
data = data.isel(olevel=0)
chl = data['CHL']
chl

anom = chl.groupby('time_counter.month') - clim
anom.name = 'anom_chl'
anom

from dask.diagnostics import ProgressBar
delayed = anom.to_netcdf('anom_chl_monthly_model.nc', compute=False)
with ProgressBar():
    delayed.compute()
