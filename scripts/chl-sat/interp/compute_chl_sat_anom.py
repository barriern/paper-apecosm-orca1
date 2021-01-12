import os.path
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from glob import glob

dirin = 'data/'

dataclim = xr.open_dataset('%s/clim_chl_monthly_obs.nc' %dirin)
clim = dataclim['clim_chl']

filelist = np.sort(glob('%s/interpolated*[0-9]*nc' %(dirin)))

for f in filelist:

    print(f)

    data = xr.open_dataset(f)
    month = data['time.month'].values[0] - 1

    tmpclim = clim.isel(month=slice(month, month + 1))
    
    count = data['count_a']

    out = data['chlor_a'] - tmpclim.values
    out = out.where(count >= (100 / 3.))
    
    out = out.rename('chl_anom')
    
    fout = os.path.basename(f)
    fout = '%s/anom_%s' %(dirin, fout)

    out.to_netcdf(fout)

