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

data = xr.open_mfdataset('%s/equatorial_ORCA1_JRAC02_CORMSK_CYC3_FINAL_FORAGE*nc' %(dirin), combine='by_coords')
time = data['time']
forage = data['FORAGE'].values   # time, dn, x, z, comm, sizes
forage[np.isnan(forage)] = 0

clim, forage = ts.get_monthly_clim(forage)
clim = np.mean(clim, axis=0).T
print(forage.shape)

forage = forage.T   # sizes, commm, z, x, dn, time
forage = sig.detrend(forage)
print(forage.shape)

dsout = xr.Dataset()
dsout['time'] = time
dsout['forage_anoms'] = (['sizes', 'comm', 'z', 'x', 'dn', 'time'], forage)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('%s/detrended_equatorial_forage_anomalies.nc' %dirout)

dsout = xr.Dataset()
dsout['forage_mean'] = (['sizes', 'comm', 'z', 'x', 'dn'], clim)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('%s/equatorial_forage_mean.nc' %dirout)
