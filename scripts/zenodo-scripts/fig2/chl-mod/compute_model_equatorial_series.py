import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os.path

mesh = xr.open_dataset('../../data/static/pacific_mesh_mask.nc').isel(z=0)
tmask = mesh['tmask'].values
lon = mesh['glamt'].values
lat = mesh['gphit'].values
e1t = mesh['e1t'].values
e2t = mesh['e2t'].values

tmask.shape, lon.shape, e1t.shape, lat.shape, e2t.shape

output = tmask.copy()
output

lonmin = 150
lonmax = -80

test = (np.abs(lat) <= 2)
test2 = (lon <= lonmax) | (lon >= lonmin)
test = (test & test2 & (tmask == 1))

ilat, ilon = np.nonzero(test == True)
output[ilat, ilon] = 2

surf = mesh['e1t'] * mesh['e2t'] * mesh['tmaskutil']
surf = surf.where(test)
surf.plot()

data = xr.open_mfdataset('../../data/nemo-pisces/*add_T*nc').isel(olevel=0)
time = data['time_counter']
data = data['CHL']
data

output = (surf * data).sum(dim=['x', 'y']) / surf.sum(dim=['x', 'y'])
output.name = 'chl'
output

from dask.diagnostics import ProgressBar
delayed = output.to_netcdf('simulated_equatorial_mean.nc', compute=False)
with ProgressBar():
    delayed.compute()
