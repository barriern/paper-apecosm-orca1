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
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

anom = False
window = 12
ystart = 1985
yend = 2016
start = '%d-01-01' %ystart
end = '%d-12-31' %yend
prefix = ''
# -

apecosm = xr.open_dataset('%sbarycenter_apecosm.nc' %prefix)
apecosm = apecosm['__xarray_dataarray_variable__']
apecosm = apecosm.sel(time=slice(start, end))
if(anom):
    apecosm = (apecosm - apecosm.mean()) / apecosm.std()
apecosm.shape

years = apecosm['time.year'].values

months = apecosm['time.month'].values

dates = years * 100 + months

ntime = len(dates)
time = np.arange(ntime)

sardara = xr.open_dataset('barycenter_sardara.nc')
sardara = sardara['catch']
sardara = sardara.sel(time=slice(ystart * 100 + 1, yend * 100 + 12))
sardara['time'].values[0]
#sardara = sardara.rolling(time=window, center=True).mean()
if(anom):
    sardara = (sardara - sardara.mean()) / sardara.std()
sardara.shape

old_sardara = xr.open_dataset('barycenter_old_sardara.nc')
old_sardara = old_sardara['catch']
old_sardara = old_sardara.sel(time=slice(ystart * 100 + 1, yend * 100 + 12))
old_sardara['time'].values[0]
#old_sardara = old_sardara.rolling(time=window, center=True).mean()
old_sardara['time']
if(anom):
    old_sardara = (old_sardara - old_sardara.mean()) / old_sardara.std()
old_sardara.shape

# +
istart = 270
istart2 = -370
toplot1 = old_sardara.values

toplot2 = sardara.values

datestr = ['%.4d-%.2d' %(y,m) for y, m in zip(years, months)]

plt.rcParams['font.size'] = 15
plt.figure(facecolor='white', figsize=(12, 8))
ax = plt.gca()
plt.plot(time, apecosm.values, label="Apecosm")
plt.grid(True)
plt.plot(time, sardara.values, label="Sardara")
plt.plot(time, old_sardara.values, label="New Sardara")
stride = 3*12
t = ax.set_xticks(time[::stride])
tl = ax.set_xticklabels(datestr[::stride], ha='right', rotation=45)
ax.set_xlim(time.min(), time.max())
plt.grid(True)
plt.legend(loc=0)
l = plt.ylabel('Longitude of biomass bary.')
plt.savefig('%sbarycenters.png' %prefix, bbox_inches='tight')
# -


