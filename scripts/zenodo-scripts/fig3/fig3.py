# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python [conda env:nbarrier2]
#     language: python
#     name: conda-env-nbarrier2-py
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
from scipy.ndimage import gaussian_filter1d
filtparams = {}
filtparams['sigma'] = 3.5
filtparams['mode'] = 'constant'
filtparams['cval'] = np.nan
filtparams['truncate'] = 4

zenodo = False
ystart = 2008
yend = 2018
start = '%d-01-01' %ystart
end = '%d-12-31' %yend
lonmin = 120
lonmax = 210
lonmax = -120 + 360
roll = False
window = 5

# +
apecosm = xr.open_dataset('apecosm/integrated_biomass_30-70cm-10N-10S.nc')
apecosm = apecosm['OOPE']
apecosm['x'] = (apecosm['x'] + 360) % 360
apecosm = apecosm.sortby(apecosm['x'])
apecosm = apecosm.sel(x=slice(lonmin, lonmax))
baryape = (apecosm['x'] * apecosm).sum(dim=['x']) / apecosm.sum(dim='x')
if(roll):
    baryape = baryape.rolling(time=window, center=True).mean()
    
apecosm = apecosm.sel(time=slice(start, end)) / (4e6 * 1e3)
# -

data1 = xr.open_dataset('sardara/regridded_catch_gear_PS_species_SKJ_1x1.nc')
data2 = xr.open_dataset('sardara/regridded_catch_gear_PS_species_YFT_1x1.nc')
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
mesh = xr.open_dataset('../data/static/pacific_mesh_mask.nc').isel(z=0, y=ilat)
lonf = mesh['glamf'].values
latf = mesh['gphif'].values
lont = mesh['glamt'].values
latt = mesh['gphit'].values

const = xr.open_dataset('../data/static/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100

compo_ape = xr.open_dataset('map_to_plot_ape.nc')
compo_ape = compo_ape['OOPE']

compo_sar = xr.open_dataset('map_to_plot_sar.nc')
compo_sar = compo_sar['catch']
compo_sar

# ## Processing the filtering/detrending of the time series

apedates = np.array([y * 100 + m for y, m in zip(baryape['time.year'], baryape['time.month'])])
iokape = np.nonzero((apedates >= 198501) & (apedates <= 201812))[0]
time = np.arange(apedates.shape[0])
apefilt = gaussian_filter1d(baryape.values.copy(), **filtparams)

sardates = barysar['time'].values
ioksar = np.nonzero((sardates >= 198501) & (sardates <= 201812))[0]
sarfilt = gaussian_filter1d(barysar.values.copy(), **filtparams)

# Extracting the filtered time-series on the common periods, and where there are no NaNs.

tempape = apefilt[iokape]
tempsar = sarfilt[ioksar]
timeape = time[iokape]
test = np.isnan(tempape) | np.isnan(tempsar)
itest = np.nonzero(test == False)[0]
tempape = tempape[itest]
tempsar = tempsar[itest]
timeape = timeape[itest]

# +
# parameters to be estimated
alpha_init = 87.0
beta_init = 1.0

# constants for optimization
x0 = np.array([alpha_init, beta_init])
bnds = ((-1., 2.), (0.8, 1.2))

def compute_series(param):
    alpha = alpha_init * param[0]
    beta = beta_init * param[1]
    shift_max = np.log(alpha + 731) * beta
    sardara_detrend = tempsar * shift_max / np.log(alpha + timeape)
    return sardara_detrend

def calc_detrend(param):
    cost = 0.
    sardara_detrend = compute_series(param)
    cost = np.sum((sardara_detrend - tempape)**2)
    return(cost)

import scipy.optimize as optimize
print('OPTIMIZATION differential evolution')
res0 = optimize.differential_evolution(calc_detrend, bounds=bnds, tol=1.E-4)
res0

print('OPTIMIZATION gradient')
res1 = optimize.minimize(calc_detrend, x0, bounds=bnds, tol=1.E-9, options={'disp': True})
res1

param0 = res0.x
param1 = res1.x

corr0 = compute_series(param0)
corr1 = compute_series(param1)

ioffset = np.nonzero(apedates == 198501)[0][0]
# -

# ## Plotting

# +
plt.rcParams['font.size'] = 15
formatter0 = LongitudeFormatter(dateline_direction_label=True)
plt.rcParams['image.cmap'] = 'RdBu_r'
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)

ymax = 1985

xloc_hov = np.arange(160, -120 + 360, 20)

plt.figure(figsize = (14, 12), facecolor='white')

############################################################ plotting barycenter time series

xxx0 = 0.02
ax = plt.axes([0.55, xxx0, 0.3, 1 - 2 * xxx0])

datebaryape = baryape['time'].values
datebaryape = np.array([d.year * 100 + d.month for d in datebaryape])

yyy = datebaryape // 100
mmm = datebaryape - 100 * yyy

iape = np.nonzero(yyy >= ymax)[0]
yyy = yyy[iape]
mmm = mmm[iape]
baryape = baryape[iape]
apefilt = apefilt[iape]

time = np.arange(len(baryape))
datestr = np.array(['%.4d-%.2d' %(y, m) for y, m in zip(yyy, mmm)])
stride = 24
ax.set_yticks(time[::stride])
ax.set_yticklabels(datestr[::stride], ha='right', rotation=45)

iline = np.nonzero(datestr=='2008-01')[0]
lw = 3

#plt.plot(np.arange(baryape.shape[0]), baryape, label='Apecosm', linewidth=0.5)
plt.plot(apefilt, np.arange(apefilt.shape[0]), label='Apecosm', linewidth=lw)
#plt.plot(np.arange(len(barysar[ioksar])), barysar[ioksar], label='Sardara', linewidth=0.5)
plt.plot(sarfilt[ioksar], np.arange(sarfilt[ioksar].shape[0]), label='Raw Obs.', linewidth=0.5 * lw, linestyle='--', color='green')
plt.plot(corr0, np.arange(corr0.shape[0]), label='Det. Obs.', linewidth=lw, color='orange')

plt.legend(ncol=1, fontsize=15, loc='upper left')
# plt.axhline(time[iline], color='k', linestyle='--')
plt.ylim(time.min(), time.max())

ax.xaxis.set_major_formatter(formatter0)
ax.set_xticks(np.arange(150, -170 + 360, 10))
plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

#ax.set_xlim(xlim)
ax.set_title('Biomass and catch barycenters')
plt.text(180, time[-10], 'd)', bbox=dicttext, ha='center', va='center')
#ax.yaxis.tick_right()
ax.grid(True, color='k', linewidth=1)

########################################################### Plotting Hovmoller diagram

pos1 = ax.get_position()

xcoord = lambda year : pos1.y0 + (year - 1985) * (pos1.height) / (2019 - 1985) 

xxx = xcoord(2008)
#ax = plt.axes([xxx, 0.43, pos1.x0 + pos1.width - xxx, 0.3])
ax = plt.axes([0.01, xxx, 0.45, pos1.y0 + pos1.height - xxx])
pos2 = ax.get_position()

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
stride = 24
ax.set_yticks(time[12::stride])
#ax.set_yticklabels(datestr[::stride], va='top', rotation=45)
plt.setp(ax.get_yticklabels(), visible=False)

ax.xaxis.set_major_formatter(formatter0)
ax.set_xticks(xloc_hov)
plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
plt.text(apecosm['x'].values[-10], time[-10],  'a)', bbox=dicttext, ha='center', va='center', zorder=3000, )

cb = plt.colorbar(cs, orientation='vertical', pad=0.05, location='left')
ax.set_xlim(140, -120+360)
xlim = ax.get_xlim()
ax.set_title('Catch (colors) \& sim. biomass (contours)')
cb.set_label('Tons (Log-scale)')
ax.grid(True, zorder=10000, color='k', linewidth=1)

############################################################ plotting catch composites

offsetx = 70
xxx0 = -0.0 + 0.03
yyy0 = 0.2 + 0.37 - 0.2
www0 = 0.45
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
gl.bottom_labels = False
gl.xlocator = mticker.FixedLocator(np.arange(-180, 180 + 20, 20))
gl.ylocator = mticker.FixedLocator(np.arange(-90, 90 + 20, 20))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

cs = ax.pcolormesh(compo_sar['lon'], compo_sar['lat'], compo_sar.values, shading='auto', transform=projin)
ccc = 3e-8
cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, color='lightgray')
ax.add_feature(cfeature.COASTLINE)
ax.set_title('NINO - NINA composites (Catches)')
ax.set_extent([130, -60 + 360 - offsetx, -20, 20], crs=projin)
plt.text(-140, 10,  'b)', bbox=dicttext, ha='center', zorder=3000, va='center', transform=projin)
#plt.text(compo_sar['lon'].values[30], compo_sar['lat'].values[10], 'Catches', bbox=dicttext, ha='center', va='center', transform=projin, zorder=3000)

#################################################################### Plot Apecosm composites

yyy0 -= 0.1 + 0.15
#xxx0 += 0.1
pos = np.array([xxx0, yyy0, www0, hhh0])

ax = plt.axes(pos, projection=ccrs.PlateCarree(central_longitude=180))
projin = ccrs.PlateCarree()

gridparams = {'crs': ccrs.PlateCarree(central_longitude=0),
              'draw_labels':True, 'linewidth':0.5,
              'color':'k', 'alpha':1, 'linestyle':'--'}
gl = ax.gridlines(**gridparams, zorder=10)
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(np.arange(-180, 180 + 20, 20))
gl.ylocator = mticker.FixedLocator(np.arange(-90, 90 + 20, 20))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

cs = ax.pcolormesh(lonf, latf, compo_ape.values[1:, 1:], transform=projin)
ax.add_feature(cfeature.LAND, color='lightgray')
ax.add_feature(cfeature.COASTLINE)
ccc = 3e-8
cs.set_clim(-ccc, ccc)
ax.set_extent([130, -60 + 360 - offsetx, -20, 20], crs=projin)
plt.text(-140, 10,  'c)', bbox=dicttext, ha='center', va='center', transform=projin, zorder=3000)
ax.set_title('NINO - NINA composites (Biomass)')

pos = np.array([xxx0, 0.03, www0, 0.02])
cax = plt.axes(pos)
cb = plt.colorbar(cs, cax=cax, orientation='horizontal', shrink=0.5)
cb.set_label('Tons/m2')

plt.savefig('gr3.jpg', bbox_inches='tight')
# -
ts1 = tempape
ts2 = corr0
np.corrcoef(ts1, ts2)[0, 1]
ts1 = (ts1 - np.mean(ts1)) / np.std(ts1)
ts2 = (ts2 - np.mean(ts2)) / np.std(ts2)
np.corrcoef(ts1, ts2)[0, 1]

# +
# import sys
# sys.path.append('../')
# from significativity import sig
# sig(ts1, ts2, use_bres=False, dof=2)

# +
# sig(ts1, ts2, use_bres=True, dof=2)
