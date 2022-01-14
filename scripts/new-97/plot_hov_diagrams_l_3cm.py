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
letters = string.ascii_letters
warnings.filterwarnings("ignore", category=RuntimeWarning) 
plt.rcParams['image.cmap'] = 'RdBu_r'
plt.rcParams['font.size'] = 15
formatter0 = LongitudeFormatter(dateline_direction_label=True)

l0 = 3
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

    time = np.arange(toplot.shape[0]) + 1

    cs = ax.pcolormesh(lon[ilon], time, toplot, shading='auto')
    cs.set_clim(cmin, cmax)
    if contour:
        cl = ax.contour(lon[ilon], time, toplot, levels=levels, colors='k', linewidths=0.5)
        cl2 = ax.contour(lon[ilon], time, toplot, levels=0, linewidths=1, colors='k')
    stride = 3
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


def compute_thetao(add_clim=True):
   
    data = xr.open_dataset('data/nino_equatorial_composites_thetao.nc')['thetao'].isel(olevel=0)
    data

    clim = xr.open_dataset('data/equatorial_full_thetao.nc')['thetao']
    clim =  clim.sel(time_counter=slice('1971-01-01', '2000-12-31')).groupby('time_counter.month').mean(dim='time_counter').isel(olevel=0)
    clim
    
    thetao = data.to_masked_array().copy()
    nyears = thetao.shape[0] // 12
    for y in range(nyears):
        index = slice(0 + y * 12, 12 + y * 12)
        thetao[index] += clim.to_masked_array()
    return thetao


# +
fig = plt.figure(figsize=(13, 14), facecolor='white')
plt.rcParams['font.size'] = 15

thetao = compute_thetao()

ilon = np.nonzero((lon >= 150) & (lon <= -90 + 360))[0]

time1 = 2
lon1 = 255

time2 = len(time) - 3
lon2 = 255

axgr = ImageGrid(fig, 111, nrows_ncols=(3, 2), axes_pad=(1.1, 0.5), cbar_pad=0.1, direction='row', aspect=False, cbar_mode="each", share_all=True)

ccc = 150
cpt = -1

cpt += 1
ax = axgr[cpt]
toplot = oope.values
toplot = toplot - toplot[0]
cs, cl = plot(ax, toplot, wstep, clim=ccc)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'B', **textprop)
cb.set_label('J/m2')
#cl = ax.contour(lon[ilon], time + 1, thetao[:, ilon], levels=[28], colors='k')
#plt.clabel(cl)

cpt += 1
toplot = (zadv + madv + zdiff + mdiff + pred + growth).cumsum(dim='time')
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep, clim=ccc)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'T', **textprop)
cb.set_label('J/m2')

cpt += 1
toplot = (growth).cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep, clim=500)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'G', **textprop)
cb.set_label('J/m2')

cpt += 1
toplot = (pred).cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep, clim=500)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'P', **textprop)
cb.set_label('J/m2')

cpt += 1
toplot = (pred + growth).cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep, clim=ccc)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'P+G', **textprop)
cb.set_label('J/m2')

cpt += 1
toplot = (zadv + madv + zdiff + mdiff).cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep, clim=ccc)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'A+D', **textprop)
cb.set_label('J/m2')

plt.savefig('hov_compo_l_%d.png' %l0, bbox_inches='tight')
# -

