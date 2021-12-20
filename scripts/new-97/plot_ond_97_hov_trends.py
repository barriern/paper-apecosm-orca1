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
#     display_name: Python [conda env:nbarrier] *
#     language: python
#     name: conda-env-nbarrier-py
# ---

# ## Imports libraries

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid
from cartopy.mpl.ticker import (LatitudeFormatter, LongitudeFormatter,
                                LatitudeLocator, LongitudeLocator)
import string
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
plt.rcParams['image.cmap'] = 'RdBu_r'
plt.rcParams['font.size'] = 15
formatter0 = LongitudeFormatter(dateline_direction_label=True)

# ## Reading weight step variable

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l', 'gpred': 'c'})
const['l'] = const['length'].values * 100
length = const['l'].values
wstep = const['weight_step']

# ## Reading OOPE hovmoller anomalies
#
# Only to check that the reconstruction of the anomalies are ok. We also shift the `x` coordinates to make them Pacific like.

data = xr.open_dataset('data/hovmoller_plot_data_oope.nc')
oope = data['OOPE']
lon = oope['x'].values
oope['x'] = (oope['x'] + 360) %360
oope = oope.sortby(oope['x'])

# ## Reconstruction of OOPE anomalies from trends
#
# First, we compute the initial OOPE values used for the integration. Since the trends are integrated from the first to  the last day of the month, we average the OOPE anomalies betweeen 1996-12 and 1997-01, in order to get close to OOPE values for 1997-01-01 (instead of 1997-01-16)

data = xr.open_dataset('data/pacific_nino97_OOPE.nc').isel(y=61).sel(time=slice('1996-12-16', '1997-01-31'))['OOPE']
dataclim = xr.open_mfdataset('data/pacific_clim_OOPE.nc').isel(y=61)['OOPE']
oope0 = data.groupby('time.month') - dataclim
oope0 = oope0.rename({'w': 'l'})
oope0['x'] = lon
oope0['x'] = (oope0['x'] + 360) %360
oope0['l'] = length
oope0 = oope0.mean(dim='time', keepdims=True)
oope0

# Now we define the list of trends to use in the integration:

varnames = [
    'growthTrend',
    'predationTrend',
    'zadv_trend',
    'zdiff_trend',
    'madv_trend',
    'mdiff_trend',
]
varnames


# Now we create the function that computes the OOPE anomalies based on trend anomalies. We first add the integrated trends to the initial OOPE. And then we compute the average betweeen step `t` and step `t+1` in order to center the OOPE in the middle of the month instead of at the end.

def compute_trend(varnames):
    for v in varnames:
        filename = 'data/hovmoller_plot_data_%s.nc' %v
        data = xr.open_dataset(filename)[v]
        if(v == varnames[0]):
            trend = data
        else:
            trend += data
    trend['x'] = (trend['x'] + 360) %360
    trend = trend.sortby(trend['x'])
    trend = trend + oope0.isel(time=0)
    
    trendnp = np.concatenate([oope0.values, trend.values], axis=0)
    trendnp = 0.5 * (trendnp[1:, ...] + trendnp[:-1, ...])
    trend.values = trendnp
    return trend
trend = compute_trend(varnames)
trend


# ## Plotting the results
#
# Now we plot the results. First, we define a function that returns the colormap limit:

def get_clim(var):
    var = var[~np.isnan(var)]
    perc = np.percentile(np.abs(var), 99)
    return -perc, perc


# Now we define the coordinates that will be used in the drawing:

lon = oope['x'].values
year = oope['time.year'].values
month = oope['time.month'].values
time = np.arange(0, oope.shape[0])
time

datestr = ['%.4d-%.2d' %(y, m) for y, m in zip(year, month)]
datestr = np.array(datestr)
datestr

letters = string.ascii_lowercase
letters


# Now we create a function that makes the drawing:

def plot(ax, toplot, wstep, contour=True, levels=None, clim=None, trend=True):
    ilon = np.nonzero((lon >= 150) & (lon <= -90 + 360))[0]
    toplot = toplot[:, ilon]
    toplot = toplot * wstep
    if clim is None:
        cmin, cmax = get_clim(toplot)
    else:
        cmin, cmax = clim
    if levels is None:
        levels = np.linspace(cmin, cmax, 11)

    time = np.arange(toplot.shape[0])

    cs = ax.pcolormesh(lon[ilon], time, toplot, shading='auto')
    cs.set_clim(cmin, cmax)
    if contour:
        cl = ax.contour(lon[ilon], time, toplot, levels=levels, colors='k', linewidths=0.5)
        cl2 = ax.contour(lon[ilon], time, toplot, levels=0, linewidths=1, colors='k')
    stride = 3
    ax.set_yticks(time[::stride])
    ax.set_yticklabels(datestr[::stride])
    ax.xaxis.set_major_formatter(formatter0)
    ax.grid(True, linewidth=0.5, color='gray', linestyle='--')
    return cs


dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)

# +
fig = plt.figure(figsize=(14, 12), facecolor='white')
plt.rcParams['font.size'] = 13

axgr = ImageGrid(fig, 111, nrows_ncols=(3, 3), axes_pad=(0.7, 0.5), cbar_pad=0.1, direction='row', aspect=False, cbar_mode="each", share_all=True)

i = -1

time1 = 2
lon1 = 255

time2 = len(time) - 3
lon2 = 255

for l in [3, 20, 90]:

    ws = float(wstep.sel(l=l, method='nearest'))
    
#     toplot = oope.sel(l=l, method='nearest').values
#     i += 1
#     cs = plot(axgr[i], toplot, ws, trend=False)
#     cb = plt.colorbar(cs, cax=axgr.cbar_axes[i])
#     cb.set_label('J/m2')
#     axgr[i].set_title('OOPE')
#     cl = cs.get_clim()

#     trend = compute_trend(varnames).sel(l=l, method='nearest')
#     toplot = trend.values
#     i += 1
#     cs = plot(axgr[i], toplot, ws)
#     cb = plt.colorbar(cs, cax=axgr.cbar_axes[i])
#     cb.set_label('J/m2')
#     if l == 3:
#         axgr[i].set_title('Total')

    trend = compute_trend(['predationTrend']).sel(l=l, method='nearest')
    toplot = trend.values
    i += 1
    cs = plot(axgr[i], toplot, ws, )
    cb = plt.colorbar(cs, cax=axgr.cbar_axes[i])
    cb.set_label('J/m2')
    if l == 3:
        axgr[i].set_title('Predation')
    axgr[i].text(lon1, time1, '%s)' %letters[i], bbox=dicttext, ha='center', va='center')
    axgr[i].text(lon2, time2, '%scm' %l, bbox=dicttext, ha='center', va='center')
        
    trend = compute_trend(['growthTrend']).sel(l=l, method='nearest')
    toplot = trend.values
    i += 1
    cs = plot(axgr[i], toplot, ws, )
    cb = plt.colorbar(cs, cax=axgr.cbar_axes[i])
    cb.set_label('J/m2')
    if(l == 3):
        axgr[i].set_title('Growth')
    axgr[i].text(lon1, time1, '%s)' %letters[i], bbox=dicttext, ha='center', va='center')
    axgr[i].text(lon2, time2, '%scm' %l, bbox=dicttext, ha='center', va='center')

    trend = compute_trend(['madv_trend', 'zadv_trend', 'mdiff_trend', 'zdiff_trend']).sel(l=l, method='nearest')
    toplot = trend.values
    i += 1
    cs = plot(axgr[i], toplot, ws)
    cb = plt.colorbar(cs, cax=axgr.cbar_axes[i])
    cb.set_label('J/m2')
    if(l == 3):
        axgr[i].set_title('Adv. + Diff')
    axgr[i].text(lon1, time1, '%s)' %letters[i], bbox=dicttext, ha='center', va='center')
    axgr[i].text(lon2, time2, '%scm' %l, bbox=dicttext, ha='center', va='center')

plt.savefig('hovmoller_anoms_oope_trends.png', bbox_inches='tight')
# -


