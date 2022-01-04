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
# -

apecosm = xr.open_dataset('barycenter_apecosm.nc')
apecosm = apecosm['__xarray_dataarray_variable__']
apecosm = apecosm.rolling(time=window, center=True).mean()
if(anom):
    apecosm = (apecosm - apecosm.mean()) / apecosm.std()

years = apecosm['time.year'].values

months = apecosm['time.month'].values

dates = years * 100 + months

ntime = len(dates)
time = np.arange(ntime)

sardara = xr.open_dataset('barycenter_sardara.nc')
sardara = sardara['catch']
sardara['time'].values[0]
sardara = sardara.rolling(time=window, center=True).mean()
if(anom):
    sardara = (sardara - sardara.mean()) / sardara.std()

old_sardara = xr.open_dataset('barycenter_old_sardara.nc')
old_sardara = old_sardara['catch']
old_sardara['time'].values[0]
old_sardara = old_sardara.rolling(time=window, center=True).mean()
old_sardara['time']
if(anom):
    old_sardara = (old_sardara - old_sardara.mean()) / old_sardara.std()

iok = np.nonzero((dates >= sardara['time'].values[0]) & (dates <= sardara['time'].values[-1]))[0]
dates[iok]

old_iok = np.nonzero((old_sardara['time'].values >= 195801) & ( old_sardara['time'].values <= 201812))[0]
old_iok

# +
istart = 270
istart2 = -370
toplot1 = old_sardara.values[old_iok][istart2:]
toplot1 = sig.detrend(toplot1) + toplot1.mean()

toplot2 = sardara.values[istart:][:-16]
toplot2 = sig.detrend(toplot2) + toplot2.mean()


plt.rcParams['font.size'] = 15
plt.figure(facecolor='white', figsize=(12, 8))
ax = plt.gca()
plt.plot(time, apecosm.values, label="Apecosm")
plt.grid(True)
plt.plot(time[iok][istart:], sardara.values[istart:], label="Sardara")
plt.plot(time[istart2:], old_sardara.values[old_iok][istart2:], label="New Sardara")
plt.plot(time[istart2:], toplot1, label="Det. New Sardara")
plt.plot(time[iok][istart:][:-16], toplot2, label="Det. Sardara")
t = ax.set_xticks(time[::stride])
tl = ax.set_xticklabels(datestr[::stride], ha='right', rotation=45)
ax.set_xlim(time.min(), time.max())
plt.grid(True)
plt.legend(loc=0)
l = plt.ylabel('Longitude of biomass bary.')
plt.savefig('barycenters.png', bbox_inches='tight')
# -

sardara.values[istart:][:-16]


