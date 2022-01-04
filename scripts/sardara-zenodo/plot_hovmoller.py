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
import matplotlib.pyplot as plt
import numpy as np
ymin = 2014
ymax = 2017

obs = xr.open_dataset('hovmoller_old_sardara.nc')
obs = obs['catch']
obs

obs = obs.sel(time=slice(ymin*100 + 1, ymax*100+12))
obs

lonobs = obs['lon']

data = xr.open_dataset('integrated_biomass_30-70cm-10N-10S.nc')
data = data['OOPE']

data = data.sel(time=slice('%d-01-01' %ymin, '%d-12-31' %ymax))
data

ntime, nx = data.shape
time = np.arange(ntime)
lon = data['x'].values
lon[lon<0] += 360

data

years = data['time.year']
years

months = data['time.month']
months

datestr = ['%.4d-%.2d' %(y, m) for y, m in zip(years, months)]
datestr

# +
iii = np.nonzero((lon >= 150) & (lon <= 260))[0]

plt.figure(facecolor='white')
ax = plt.gca()

toplot = data.values.copy()
toplot = np.log10(toplot, where=toplot>0, out=toplot)
toplot = np.ma.masked_where(toplot==0, toplot)
cl = plt.contour(lon[iii], time, toplot[:, iii], 15, colors='k', linewidths=1)
plt.clabel(cl)

toplot = obs.values.copy()
toplot = np.ma.masked_where(toplot==0, toplot)
toplot = np.log10(toplot, where=toplot>0, out=toplot)
cs = plt.pcolormesh(lonobs, time, toplot, shading='auto')
plt.colorbar(cs)

stride = 6
ax.set_yticks(time[::stride])
ax.set_yticklabels(datestr[::stride], va='top', rotation=45)


plt.gca().set_xlim(150, 260)
plt.savefig('hovmoller_bis.png', bbox_inches='tight')
# -


