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
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker

zenodo = False
ystart = 2008
yend = 2018
start = '%d-01-01' %ystart
end = '%d-12-31' %yend
lonmin = 120
lonmax = 210
lonmax = -120 + 360
roll = True
window = 5

# +
apecosm = xr.open_dataset('integrated_biomass_30-70cm-10N-10S.nc')
apecosm = apecosm['OOPE']
apecosm['x'] = (apecosm['x'] + 360) % 360
apecosm = apecosm.sortby(apecosm['x'])
apecosm = apecosm.sel(x=slice(lonmin, lonmax))
baryape = (apecosm['x'] * apecosm).sum(dim=['x']) / apecosm.sum(dim='x')
if(roll):
    baryape = baryape.rolling(time=window, center=True).mean()
    
apecosm = apecosm.sel(time=slice(start, end)) / 4e6
# -

if zenodo:
    catch = xr.open_dataset('regridded_catch_gear_PS.nc')
    catch = catch.where(abs(catch['lat']) <= 10, drop=True)
    catch = catch['catch'].sum(dim=['species', 'lat'])
    catch = catch.sel(lon=slice(lonmin, lonmax))
    catch = catch.where(catch != 0)
    barysar = (catch['lon'] * catch).sum(dim=['lon']) / catch.sum(dim='lon')
else:
    data1 = xr.open_dataset('../sardara/data/regridded_catch_gear_PS_species_SKJ_1x1.nc')
    data2 = xr.open_dataset('../sardara/data/regridded_catch_gear_PS_species_YFT_1x1.nc')
    catch = data1['catch'] + data2['catch']
    catch = catch.where(abs(catch['lat']) <= 10, drop=True)
    catch = catch.where(catch != 0)
    catch = catch.sum(dim='lat')
    catch = catch.sel(lon=slice(lonmin, lonmax))
    barysar = (catch['lon'] * catch).sum(dim=['lon']) / catch.sum(dim='lon')
if(roll):
    barysar = barysar.rolling(time=window, center=True).mean()
catch = catch.sel(time=slice(ystart * 100, yend * 100 + 12))

ilat = slice(None, -3)
mesh = xr.open_dataset('../new-97/data/pacific_mesh_mask.nc').isel(z=0, y=ilat)
lonf = mesh['glamf'].values
latf = mesh['gphif'].values
lont = mesh['glamt'].values
latt = mesh['gphit'].values

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100

compo_ape = xr.open_dataset('map_to_plot_ape.nc')
compo_ape = compo_ape['OOPE']

compo_sar = xr.open_dataset('map_to_plot_sar.nc')
compo_sar = compo_sar['catch']
compo_sar

# +
plt.rcParams['font.size'] = 15
formatter0 = LongitudeFormatter(dateline_direction_label=True)
plt.rcParams['image.cmap'] = 'RdBu_r'

plt.figure(figsize = (12, 14), facecolor='white')

ax = plt.axes([0.05, 0.5, 0.3, 0.4])

apeyears = apecosm['time.year'].values
apemonths = apecosm['time.month'].values
apedates = apeyears * 100 + apemonths
time = np.arange(len(apedates))

toplot = np.log10(apecosm.values)
cl = plt.contour(apecosm['x'].values, time, toplot, 6, linewidths=1, colors='k')
plt.clabel(cl)

toplot = np.log10(catch.values, where=catch>0, out=np.zeros(catch.values.shape))
toplot = np.ma.masked_where(toplot == 0, toplot)
cs = plt.pcolormesh(catch['lon'].values, time, toplot, shading='auto', cmap=plt.cm.Spectral_r)
cs.set_clim(0, 4.5)

datestr = ['%.4d-%.2d' %(y, m) for y, m in zip(apeyears, apemonths)]
stride = 12
ax.set_yticks(time[::stride])
ax.set_yticklabels(datestr[::stride], va='top', rotation=45)

ax.xaxis.set_major_formatter(formatter0)
ax.set_xticks(np.arange(160, -120 + 360, 20))
ax.set_xticks(np.arange(140, -120 + 360, 20))
plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

cb = plt.colorbar(cs)
ax.set_xlim(140, -120+360)
xlim = ax.get_xlim()

############################################################ plotting lower panel

pos1 = ax.get_position()
offset = 0.45
pos2 = [pos1.x0, pos1.y0 -offset,  pos1.width, pos1.height]
pos2 = [pos1.x0 + offset, pos1.y0,  pos1.width, pos1.height]

ax = plt.axes(pos2)

barysartemp = barysar.sel(time=slice(198501, 201812))
datebaryape = baryape['time'].values
datebaryape = [d.year * 100 + d.month for d in datebaryape]
datebarysar = barysartemp['time'].values

yyy = datebarysar // 100
mmm = datebarysar - 100 * yyy

test = (datebaryape >= datebarysar[0]) & (datebaryape <= datebarysar[-1])
iok = np.nonzero(test)[0]

time = np.arange(len(datebarysar))
datestr = np.array(['%.4d-%.2d' %(y, m) for y, m in zip(yyy, mmm)])
stride = 4*12
ax.set_yticks(time[::stride])
ax.set_yticklabels(datestr[::stride], va='top', rotation=45)

iline = np.nonzero(datestr=='2008-01')[0]

plt.plot(baryape[iok], time, label='Apecosm')
plt.plot(barysartemp, time, label='Sardara')
plt.legend()
plt.axhline(time[iline], color='k', linestyle='--')
plt.ylim(time.min(), time.max())

ax.xaxis.set_major_formatter(formatter0)
ax.set_xticks(np.arange(160, -120 + 360, 20))
plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

cb.set_label('Log(MT)')
ax.set_xlim(xlim)
ax.set_title('Catches and Biomass')

ax.set_title('Biomass and catch barycenters')

############################################################# plotting right panel

xxx0 = -0.05
yyy0 = 0.0
www0 = 0.9
hhh0 = 0.2
pos = np.array([xxx0, yyy0, www0, hhh0])

ax = plt.axes(pos, projection=ccrs.PlateCarree(central_longitude=180))
projin = ccrs.PlateCarree()

gridparams = {'crs': ccrs.PlateCarree(central_longitude=0),
              'draw_labels':True, 'linewidth':0.5,
              'color':'k', 'alpha':1, 'linestyle':'--'}
gl = ax.gridlines(**gridparams, zorder=10)
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(np.arange(-180, 180 + 40, 40))
gl.ylocator = mticker.FixedLocator(np.arange(-90, 90 + 20, 20))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

cs = ax.pcolormesh(compo_sar['lon'], compo_sar['lat'], compo_sar.values, shading='auto', transform=projin)
space = 50
#cs = ax.pcolormesh(lonf, latf, toplot2[1:, 1:], linewidths=1, transform=projin)
#cl = plt.tricontour(lonout, latout, toplot1d, colors='k', linewidths=0.5, levels=np.arange(-400 - space, 400 + space, space))
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
cb = plt.colorbar(cs, shrink=0.7, location='right', pad=0.07)
cb.set_label('Catch anoms. (MT/m2)')
ax.set_title('OND15 - 0.5 x (OND12 + OND13)')
ax.set_extent([130, -60 + 360, -40, 40], crs=projin)

#################################################################### Apecosm

yyy0 += hhh0 + 0.03
pos = np.array([xxx0, yyy0, www0, hhh0])

ax = plt.axes(pos, projection=ccrs.PlateCarree(central_longitude=180))
projin = ccrs.PlateCarree()

gridparams = {'crs': ccrs.PlateCarree(central_longitude=0),
              'draw_labels':True, 'linewidth':0.5,
              'color':'k', 'alpha':1, 'linestyle':'--'}
gl = ax.gridlines(**gridparams, zorder=10)
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(np.arange(-180, 180 + 40, 40))
gl.ylocator = mticker.FixedLocator(np.arange(-90, 90 + 20, 20))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

cs = ax.pcolormesh(lonf, latf, compo_ape.values[1:, 1:], transform=projin)
#space = 50
#cs = ax.pcolormesh(lonf, latf, toplot2[1:, 1:], linewidths=1, transform=projin)
#cl = plt.tricontour(lonout, latout, toplot1d, colors='k', linewidths=0.5, levels=np.arange(-400 - space, 400 + space, space))
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
cb = plt.colorbar(cs, shrink=0.7, location='right', pad=0.07)
ccc = 3e-8
cs.set_clim(-ccc, ccc)
cb.add_lines(cl)
cb.set_label('Biomass anoms. (MT/m2)')
#ax.set_title('OND15 - 0.5 x (OND12 + OND13)')
ax.set_extent([130, -60 + 360, -40, 40], crs=projin)

plt.savefig('plot_validation_apecosm.png', bbox_inches='tight')
# -

