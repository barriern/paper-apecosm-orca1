import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os.path

data = xr.open_dataset('../data/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M-m_1641569308380.nc')
data = data['sla']
data

latmin = -5
latmax = 5
lonmax = -120 + 360
lonmin = -170 + 360

testlon = (data['longitude'] <= lonmax) & (data['longitude'] >= lonmin)
data = data.where(testlon, drop=True)
testlat = (abs(data['latitude']) <= latmax)
data = data.where(testlat, drop=True)
data.isel(time=0).plot()

weights = np.cos(np.deg2rad(data['latitude']))
weights
weights, data = xr.broadcast(weights, data)

output = (weights * data).sum(dim=['latitude', 'longitude']) / weights.sum(dim=['latitude', 'longitude'])
output.name = 'ssh'

from dask.diagnostics import ProgressBar
delayed = output.to_netcdf('../data/ssh_obs_equatorial_mean.nc', compute=False)
with ProgressBar():
    delayed.compute()
