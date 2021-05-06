import xarray as xr
import numpy as np
import os.path 
from datetime import datetime
from glob import glob
import cftime
from datetime import datetime
import apecosm.ts as ts
import scipy.signal as sig

units = 'seconds since 1950-01-01 00:00:00'
calendar='noleap'

latmax = 2

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data/'
dirout = dirin

for prefix in ['add_T']:

    data = xr.open_mfdataset('%s/equatorial_*%s*nc' %(dirin, prefix), combine='by_coords')
    time = data['time_counter']

    var = data.variables
    varout = [v for v in var if ((data[v].ndim == 3) & ('e3' not in v))]

    for v in varout:
        
        print(":::::::::::::::::::::::::: processing variable ", v)
        forage = data[v].values   # time, z, x
        forage[np.isnan(forage)] = 0

        clim, forage = ts.get_monthly_clim(forage)
        clim = np.mean(clim, axis=0).T  # x, z

        forage = forage.T   # x, z, time
        forage = sig.detrend(forage)
        print(forage.shape)

        dsout = xr.Dataset()
        dsout['time_counter'] = time
        dsout['%s_anoms' %v] = (['x', 'z', 'time_counter'], forage)
        dsout.attrs['file'] = os.path.realpath(__file__)
        dsout.to_netcdf('%s/detrended_equatorial_%s_anomalies.nc' %(dirout, v))

        dsout = xr.Dataset()
        dsout['%s_mean' %v] = (['x', 'z'], clim)
        dsout.attrs['file'] = os.path.realpath(__file__)
        dsout.to_netcdf('%s/equatorial_%s_mean.nc' %(dirout, v))
