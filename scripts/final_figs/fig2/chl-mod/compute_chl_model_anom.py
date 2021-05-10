import os.path
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from glob import glob

dirin = 'data/'

datac = xr.open_dataset('%s/clim_chl_monthly_model.nc' %dirin)
clim = datac['clim_chl'].to_masked_array()

data = xr.open_mfdataset("data/CHL*nc", combine="by_coords")
data = data.isel(olevel=0)
chl = data['NCHL'] + data['DCHL']

ntime, nlat, nlon = chl.shape
nyears = ntime // 12

index = np.arange(12)

anom = np.zeros((ntime, nlat, nlon), dtype=np.float)

for p in range(nyears):

    anom[index, :, :] = chl.isel(time_counter=index).to_masked_array() - clim
    index += 12

dsout = xr.Dataset()
dsout['anom_chl'] = (['time', 'y', 'x'], anom)
dsout['time'] = data['time_counter']

dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('%s/anom_chl_monthly_model.nc' %dirin)
