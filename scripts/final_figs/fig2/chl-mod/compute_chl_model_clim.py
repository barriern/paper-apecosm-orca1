import os.path
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from glob import glob

dirin = 'data/'

data = xr.open_mfdataset("/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/JRA_CO2/*add_T*nc", combine="by_coords")
data = data.isel(olevel=0)
chl = data['NCHL'] + data['DCHL']
chl = chl.sel(time_counter=slice('1998-01-01', '2018-12-31'))
chl

clim = chl.groupby('time_counter.month').mean(dim='time_counter')
clim.name = 'clim_chl'

# +
delayed = clim.to_netcdf('data/clim_chl_monthly_model.nc', compute=False)

from dask.diagnostics import ProgressBar
with ProgressBar():
    delayed.compute()
# -


