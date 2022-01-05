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
from cartopy.mpl.ticker import (LatitudeFormatter, LongitudeFormatter,
                                LatitudeLocator, LongitudeLocator)

zenodo = False
ystart = 1990
yend = 2018
start = '%d-01-01' %ystart
end = '%d-12-31' %yend
lonmin = 120
lonmax = 210
lonmax = -120 + 360
roll = False
# -

apecosm = xr.open_dataset('integrated_biomass_30-70cm-10N-10S.nc')
apecosm = apecosm['OOPE']
apecosm['x'] = (apecosm['x'] + 360) % 360
apecosm = apecosm.sortby(apecosm['x'])
apecosm = apecosm.sel(x=slice(lonmin, lonmax))
baryape = (apecosm['x'] * apecosm).sum(dim=['x']) / apecosm.sum(dim='x')
if(roll):
    baryape = baryape.rolling(time=12, center=True).mean()
baryape = baryape.sel(time=slice(start, end))
apecosm = apecosm.sel(time=slice(start, end)) / 4e6

apeyears = apecosm['time.year'].values
apemonths = apecosm['time.month'].values
apedates = apeyears * 100 + apemonths
time = np.arange(len(apedates))
apedates

# +
if zenodo:
    catch = xr.open_dataset('regridded_catch_gear_PS.nc')
    catch = catch.where(abs(catch['lat']) <= 10, drop=True)
    catch = catch['catch'].sum(dim=['species', 'lat'])
    catch = catch.sel(lon=slice(lonmin, lonmax))
    barysar = (catch['lon'] * catch).sum(dim=['lon']) / catch.sum(dim='lon')
else:
    data1 = xr.open_dataset('../sardara/data/regridded_catch_gear_PS_species_SKJ_1x1.nc')
    data2 = xr.open_dataset('../sardara/data/regridded_catch_gear_PS_species_YFT_1x1.nc')
    catch = data1['catch'] + data2['catch']
    catch = catch.where(abs(catch['lat']) <= 10, drop=True)
    catch = catch.sum(dim='lat')
    catch = catch.sel(lon=slice(lonmin, lonmax))
    barysar = (catch['lon'] * catch).sum(dim=['lon']) / catch.sum(dim='lon')
if(roll):
    barysar = barysar.rolling(time=12, center=True).mean()

catch = catch.sel(time=slice(ystart * 100, yend * 100 + 12))
barysar = barysar.sel(time=slice(ystart * 100, yend * 100 + 12))

# +
plt.rcParams['font.size'] = 15
formatter0 = LongitudeFormatter(dateline_direction_label=True)

plt.figure(figsize = (8, 12), facecolor='white')
ax = plt.gca()

toplot = np.log10(apecosm.values)
cl = plt.contour(apecosm['x'].values, time, toplot, 6, linewidths=1, colors='k')
plt.plot(baryape, time, color='green')
plt.clabel(cl)

toplot = np.log10(catch.values, where=catch>0, out=np.zeros(catch.values.shape))
toplot = np.ma.masked_where(toplot == 0, toplot)
cs = plt.pcolormesh(catch['lon'].values, time, toplot, shading='auto', cmap=plt.cm.Spectral_r)
cs.set_clim(0, 4.5)

plt.plot(barysar, time, color='steelblue', label='Sardara', linewidth=3)
plt.plot(baryape, time, color='k', label='Apecosm', linewidth=3)

plt.legend(loc=0)

datestr = ['%.4d-%.2d' %(y, m) for y, m in zip(apeyears, apemonths)]
stride = 24
ax.set_yticks(time[::stride])
ax.set_yticklabels(datestr[::stride], va='top', rotation=45)

ax.xaxis.set_major_formatter(formatter0)
ax.set_xticks(np.arange(140, -120 + 360, 20))

cb = plt.colorbar(cs)
cb.set_label('Catches (Log(MT))')
plt.savefig('plot_validation_apecosm.png', bbox_inches='tight')
# -

test1 = baryape.values.copy()
test2 = barysar.values.copy()
temp = (np.isnan(test1) | np.isnan(test2))
np.corrcoef(test1[~temp], test2[~temp])

barysar


