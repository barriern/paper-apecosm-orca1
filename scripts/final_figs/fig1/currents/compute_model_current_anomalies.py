# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Python [conda env:nbarrier] *
#     language: python
#     name: conda-env-nbarrier-py
# ---

# +
import xarray as xr
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

pattern = 'speed_V'
varname = 'vo'

dirin = 'data/'
mesh = xr.open_dataset('%s/pacific_mesh_mask.nc' %dirin).isel(z=0)
mesh
# -

surf = mesh['e1t'] * mesh['e2t'] * mesh['tmaskutil']
surf.plot()

lon = mesh['glamt']
lat = mesh['gphit']
lon = (lon + 360) % 360
nino34 = (abs(lat) <= 5) & (lon >= 190) & (lon <= 240) & (mesh['tmaskutil'] == 1)
nino34.plot()

dirin = '/home1/scratch/nbarrier'
data = xr.open_mfdataset('%s/*%s*nc' %(dirin, pattern)).isel(olevel=0)
data = data[varname]
data

ts = (data.where(nino34) * surf).sum(dim=['x', 'y']) / surf.sum(dim=['x', 'y'])
delayed = ts.to_netcdf('data/model_nino_34_%s.nc' %varname, compute=False)
with ProgressBar():
    delayed.compute()

data_clim = data.groupby('time_counter.month').mean(dim='time_counter')
data_anom = data.groupby('time_counter.month') - data_clim
delayed = data.to_netcdf('data/model_anoms_%s.nc' %varname, compute=False)
with ProgressBar():
    delayed.compute()
