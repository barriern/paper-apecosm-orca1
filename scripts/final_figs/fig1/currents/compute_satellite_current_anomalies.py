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

varname = 'vo'

data = xr.open_dataset('data/dataset-%s-rep-monthly.nc' %varname).isel(depth=0)
data = data[varname].chunk({'longitude' : 400})
data
# -

lon = data['longitude']
lon

lon = (lon + 360) % 360
lon

data['longitude'] = lon
data = data.sortby(data['longitude'])
data.isel(time=0).plot()

mask = (data['longitude'] >= 120) & (data['longitude'] <= 300)
mask

output = data.where(mask == True, drop=True)
output.isel(time=0).plot()

nino34 = (abs(output['latitude']) <= 5) & (output['longitude'] >= 190) & (output['longitude'] <= 240)
ts = output.where(nino34).mean(dim=['longitude', 'latitude'])
delayed = ts.to_netcdf('satellite_nino_34_%s.nc' %varname, compute=False)
with ProgressBar():
    delayed.compute()

output_clim = output.groupby('time.month').mean(dim='time')
output_anom = output.groupby('time.month') - output_clim
delayed = output.to_netcdf('satellite_anoms_%s.nc' %varname, compute=False)
with ProgressBar():
    delayed.compute()


