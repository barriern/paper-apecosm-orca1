import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import apecosm.ts as ts
import os.path

mesh = xr.open_dataset('../../../data/mesh_mask_eORCA1_v2.2.nc').isel(t=0)
tmask = mesh['tmask'].values[0]
lon = mesh['glamt'].values
lat = mesh['gphit'].values
e1t = mesh['e1t'].values
e2t = mesh['e2t'].values

output = tmask.copy()
output

latmin = -5
latmax = 5
lonmax = -120
lonmin = -170

test = (lat <= latmax) & (lat >= latmin)
test = test & (lon<=lonmax) & (lon>=lonmin)
test = test & (tmask == 1)
xrsurf = mesh['e1t'] * mesh['e2t'] * mesh['tmaskutil']
xrsurf = xrsurf.where(test)
xrsurf.plot()

ilat, ilon = np.nonzero(test == True)
output[ilat, ilon] = 2

surf = mesh['e1t'] * mesh['e2t'] * mesh['tmaskutil']
surf = surf.where(test)
surf.plot()

data = xr.open_mfdataset('/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/JRA_CO2/*ssh_T*nc')
time = data['time_counter']
data = data['zos']
data

output = (surf * data).sum(dim=['x', 'y']) / surf.sum(dim=['x', 'y'])
output.name = 'ssh'
output

from dask.diagnostics import ProgressBar
delayed = output.to_netcdf('../data/ssh_simulated_equatorial_mean.nc', compute=False)
with ProgressBar():
    delayed.compute()


