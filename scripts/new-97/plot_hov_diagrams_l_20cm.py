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

# +
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

l0 = 20
# -

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100
const

wstep = const['weight_step'].sel(l=l0, method='nearest')
wstep = float(wstep)
wstep


# +
def read_variable(varname):
    
    data = xr.open_dataset('data/nino_equatorial_composites_%s.nc' %varname)[varname]
    data['l'] = const['l']
    data = data.sel(l=l0, method='nearest')
    return data

oope = read_variable('OOPE')
oope
# -

lon = oope['x'].values
time = oope['time'].values

pred = read_variable('predationTrend')
growth = read_variable('growthTrend')
zadv = read_variable('zadv_trend')
madv = read_variable('madv_trend')
zdiff = read_variable('zdiff_trend')
mdiff = read_variable('mdiff_trend')


def get_clim(var):
    var = var[~np.isnan(var)]
    perc = np.percentile(np.abs(var), 99)
    return -perc, perc


def plot(ax, toplot, wstep, contour=True, levels=None, clim=None, trend=True):
    ilon = np.nonzero((lon >= 150) & (lon <= -90 + 360))[0]
    toplot = toplot[:, ilon]
    toplot = toplot * wstep
    if clim is None:
        cmin, cmax = get_clim(toplot)
    else:
        cmin, cmax = [-clim, clim]
    if levels is None:
        levels = np.linspace(cmin, cmax, 11)

    time = np.arange(toplot.shape[0])

    cs = ax.pcolormesh(lon[ilon], time, toplot, shading='auto')
    cs.set_clim(cmin, cmax)
    if contour:
        cl = ax.contour(lon[ilon], time, toplot, levels=levels, colors='k', linewidths=0.5)
        cl2 = ax.contour(lon[ilon], time, toplot, levels=0, linewidths=1, colors='k')
    stride = 3
    #ax.set_yticks(time[::stride])
    #ax.set_yticklabels(datestr[::stride], rotation=45, va='top')
    ax.xaxis.set_major_formatter(formatter0)
    ax.grid(True, linewidth=0.5, color='gray', linestyle='--')
    labels = ['180', '-150', '-120', '-90']
    xticks = np.array([float(l) for l in labels])
    xticks[xticks < 0] += 360
    ax.set_xticks(xticks)
    plt.setp(ax.get_xticklabels(), ha='right', rotation=45)
    ax.set_ylim(time.min(), time.max())
    ax.set_ylabel('Month')
    return cs, None


dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
textprop = {}
textprop['bbox'] = dicttext
textprop['ha'] = 'center'
textprop['va'] = 'center'

# +
fig = plt.figure(figsize=(13, 14), facecolor='white')
plt.rcParams['font.size'] = 13

time1 = 2
lon1 = 255

time2 = len(time) - 3
lon2 = 255

axgr = ImageGrid(fig, 111, nrows_ncols=(3, 2), axes_pad=(0.9, 0.5), cbar_pad=0.1, direction='row', aspect=False, cbar_mode="each", share_all=True)

cpt = 0
ax = axgr[cpt]
toplot = oope.values
toplot = toplot - toplot[0]
cs, cl = plot(ax, toplot, wstep, clim=15)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'B', **textprop)
#ax.text(lon1, time1, 'a)', **textprop)

cpt = 1
toplot = (growth + pred + zadv + madv + zdiff + mdiff).cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep, clim=15)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'T', **textprop)

cpt = 2
toplot = (growth).cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep, clim=15)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'G', **textprop)

cpt = 3
toplot = (pred).cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep, clim=15)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'P', **textprop)

cpt = 4
toplot = (pred + growth).cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep, clim=15)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'P+G', **textprop)

cpt = 5
toplot = (zadv + madv + zdiff + mdiff).cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep, clim=15)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'A+D', **textprop)
plt.savefig('hov_compo_l_%d.png' %l0)
# -

