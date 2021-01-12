import os.path
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from glob import glob

dirin = 'data/'

data = xr.open_mfdataset("data/CHL*nc", combine="by_coords")
#data = xr.open_mfdataset("data/CHL*195*nc", combine="by_coords")
data = data.isel(olevel=0)
chl = data['NCHL'] + data['DCHL']

nmonths = 12

ntime, nlat, nlon = chl.shape
nyears = ntime // 12

index = np.arange(12)

clim = np.zeros((nmonths, nlat, nlon), dtype=np.float)

for p in range(nyears):
    print(p, nyears, index)

    clim += chl.isel(time_counter=index).to_masked_array()
    index += 12

clim /= nyears

dsout = xr.Dataset()
dsout['clim_chl'] = (['month', 'y', 'x'], clim)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('%s/clim_chl_monthly_model.nc' %dirin)
