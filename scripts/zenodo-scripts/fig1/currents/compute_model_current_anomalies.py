# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python [conda env:nbarrier]
#     language: python
#     name: conda-env-nbarrier-py
# ---

# +
import xarray as xr
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os

varname = 'vo'
if varname == 'uo':
    grid = 'speed_U'
    maskvar = 'umaskutil'
else:
    grid = 'speed_V'
    maskvar = 'tmaskutil'

dirin = os.path.join('..', '..', 'data')
mesh = xr.open_dataset(os.path.join(dirin, 'static', 'pacific_mesh_mask.nc')).isel(z=0)
mesh
# -

maskvar, grid

surf = mesh['e1t'] * mesh['e2t'] * mesh[maskvar]
surf.plot()

lon = mesh['glamt']
lat = mesh['gphit']
lon = (lon + 360) % 360
nino34 = (abs(lat) <= 5) & (lon >= 190) & (lon <= 240) & (mesh['tmaskutil'] == 1)
nino34.plot()

data = xr.open_mfdataset(os.path.join(dirin, 'nemo-pisces', '*%s*' %grid)).isel(olevel=0)
data = data[varname]
data

ts = (data.where(nino34) * surf).sum(dim=['x', 'y']) / surf.where(nino34).sum(dim=['x', 'y'])
delayed = ts.to_netcdf('model_nino_34_%s.nc' %varname, compute=False)
with ProgressBar():
    delayed.compute()

data_clim = data.groupby('time_counter.month').mean(dim='time_counter')
data_anom = data.groupby('time_counter.month') - data_clim
delayed = data.to_netcdf('model_anoms_%s.nc' %varname, compute=False)
with ProgressBar():
    delayed.compute()


