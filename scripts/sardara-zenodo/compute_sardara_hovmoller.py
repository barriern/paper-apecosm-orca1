# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.pyplot as plt

data = xr.open_dataset('regridded_catch_gear_PS.nc')
data

lon = data['lon'].values
lat = data['lat'].values
time = data['time'].values
year = time // 100
month = time - 100*year

catch = data['catch'].sum(dim='species')

date = data['time'].values
date[:24]

iok = np.nonzero((date >= 199701) & (date <=199812))[0]
iok

subcatch = catch.isel(time=iok)
subcatch

subcatch = subcatch.where(abs(subcatch['lat']) <= 10)
subcatch

ntime = subcatch.shape[0]
time = np.arange(ntime)
time

lon = subcatch['lon'].values
lon

output = subcatch.sum(dim='lat')
toplot = output.values
toplot = np.log10(toplot, where=toplot >0, out=toplot)
toplot = np.ma.masked_where(toplot==0, toplot)
cs = plt.pcolormesh(lon, time, toplot, shading='auto')
plt.colorbar(cs)
plt.gca().set_xlim(150, 260)

toplot

output

output.to_netcdf('hovmoller_sardara.nc')

output


