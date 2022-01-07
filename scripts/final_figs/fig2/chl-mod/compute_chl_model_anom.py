import os.path
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from glob import glob

dirin = 'data/'

datac = xr.open_dataset('%s/clim_chl_monthly_model.nc' %dirin)
clim = datac['clim_chl']
clim

data = xr.open_mfdataset("/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/JRA_CO2/*add_T*nc", combine="by_coords")
data = data.isel(olevel=0)
chl = data['NCHL'] + data['DCHL']
chl

anom = chl.groupby('time_counter.month') - clim
anom.name = 'anom_chl'
anom

from dask.diagnostics import ProgressBar
delayed = anom.to_netcdf('%s/anom_chl_monthly_model.nc' %dirin, compute=False)
with ProgressBar():
    delayed.compute()


