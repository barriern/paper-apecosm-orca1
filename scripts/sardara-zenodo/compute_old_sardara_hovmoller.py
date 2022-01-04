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

import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.pyplot as plt

data1 = xr.open_dataset('../sardara/data/regridded_catch_gear_PS_species_SKJ_1x1.nc')
data1
data2 = xr.open_dataset('../sardara/data/regridded_catch_gear_PS_species_YFT_1x1.nc')
data2 
catch = data1['catch'] + data2['catch']

lon = catch['lon'].values
lat = catch['lat'].values
time = catch['time'].values
year = time // 100
month = time - 100*year

dates = catch['time'].values
dates[:24]

catch = catch.where(abs(catch['lat']) <= 10)
catch

ntime = catch.shape[0]
time = np.arange(ntime)
time

lon = catch['lon'].values
lon

output = catch.sum(dim='lat')
output

output.to_netcdf('hovmoller_old_sardara.nc')
