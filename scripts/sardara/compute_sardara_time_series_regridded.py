# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import cftime
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

species = 'SKJ'
gear = 'PS'
res = 1
latmax = 2

data = xr.open_dataset('data/regridded_catch_gear_%s_species_%s_%dx%d.nc' %(gear, species, res, res), decode_times=False)
data = data['catch']
# -

lon = data['lon']
lat = data['lat']
lon.values

datawest = data.where((abs(lat) < latmax) & (lon <= 180))
datawest.isel(time=-1).plot()

dataeast = data.where((abs(lat) < latmax) & (lon > 180))
dataeast.isel(time=-1).plot()

tseast = dataeast.sum(dim=['lon', 'lat'])
tseast

tswest = datawest.sum(dim=['lon', 'lat'])

tswest.plot()

tseast.plot()

dsout = xr.Dataset()
dsout['time'] = data['time']
dsout

dsout['tseast'] = tseast

dsout['tswest'] = tswest

dsout.attrs['latmax'] = latmax

dsout.to_netcdf('sardara_time_series_gear_%s_species_%s_latmax_%d.nc' %(species, gear, latmax))


