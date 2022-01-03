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
import matplotlib.pyplot as plt
import numpy as np

obs = xr.open_dataset('hovmoller_sardara.nc')
obs = obs['catch']
obs

lonobs = obs['lon']

data = xr.open_dataset('hovmoller_apecosm.nc')
data = data['__xarray_dataarray_variable__']
data

data = data.sel(time=slice('1997-01-01', '1998-12-31'))
data

ntime, nx = data.shape
time = np.arange(ntime)
lon = data['x'].values
lon[lon<0] += 360

data

# +
iii = np.nonzero((lon >= 150) & (lon <= 260))[0]

toplot = data.values.copy()
toplot = np.log10(toplot, where=toplot >0, out=toplot)
toplot = np.ma.masked_where(toplot==0, toplot)
cl = plt.contour(lon[iii], time, toplot[:, iii], 15, colors='k', linewidths=1)
#cs = plt.pcolormesh(lon[iii], time, toplot[:, iii], shading='auto')
#plt.colorbar(cs)
plt.clabel(cl)

toplot = obs.values.copy()
# toplot = np.log10(toplot, where=toplot>0, out=toplot)
toplot = np.ma.masked_where(toplot==0, toplot)
cs = plt.pcolormesh(lonobs, time, toplot, shading='auto')
plt.colorbar(cs)

plt.gca().set_xlim(150, 260)
# -

data.values


