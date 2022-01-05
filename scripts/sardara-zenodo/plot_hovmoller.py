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
import matplotlib.pyplot as plt
import numpy as np
from cartopy.mpl.ticker import (LatitudeFormatter, LongitudeFormatter,
                                LatitudeLocator, LongitudeLocator)

ymin = 2014
ymax = 2017
prefix = 'school_'
# -

obs = xr.open_dataset('hovmoller_old_sardara.nc')
obs = obs['catch']
obs

obs = obs.sel(time=slice(ymin*100 + 1, ymax*100+12))
obs

lonobs = obs['lon']

data = xr.open_dataset('%sintegrated_biomass_30-70cm-10N-10S.nc' %prefix)
data = data['OOPE'] / 4.e6

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
formatter0 = LongitudeFormatter(dateline_direction_label=True)
plt.rcParams['font.size'] = 15

plt.figure(facecolor='white', figsize=(12, 8))
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
cb = plt.colorbar(cs)
cb.set_label('Sardara catchs (Log(Metric Tons))')

stride = 3
ax.set_yticks(time[::stride])
ax.set_yticklabels(datestr[::stride], va='top', rotation=45)
ax.xaxis.set_major_formatter(formatter0)

plt.gca().set_xlim(150, 260)
plt.savefig('%shovmoller_bis.png' %prefix, bbox_inches='tight')
