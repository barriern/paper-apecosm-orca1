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

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid
from cartopy.mpl.ticker import (LatitudeFormatter, LongitudeFormatter,
                                LatitudeLocator, LongitudeLocator)
import string
import warnings
letters = string.ascii_letters
letters = [l + ')' for l in letters]
warnings.filterwarnings("ignore", category=RuntimeWarning) 
plt.rcParams['image.cmap'] = 'RdBu_r'
plt.rcParams['font.size'] = 15
formatter0 = LongitudeFormatter(dateline_direction_label=True)

mesh = xr.open_dataset('equatorial_mesh_mask.nc')
mesh = mesh.rename({'z': 'olevel'})
lon = mesh['x'].values
ilon = np.nonzero((lon >= 150) & (lon <= -90 + 360))[0]
mesh = mesh.isel(x=ilon)
lon = mesh['x'].values
depth = mesh['gdept_1d'].isel(x=0)
depth
mesh['olevel'] = depth
mesh

const = xr.open_dataset('../data/static/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100
const

wstep = const['weight_step']
wstep = wstep.sel(l=[3, 20, 90], method='nearest')


# +
def read_variable(varname):
    
    data = xr.open_dataset('../fig4/nino_equatorial_composites_%s.nc' %varname)[varname].isel(x=ilon)
    return data

oope = read_variable('OOPE')
oope['l'] = const['l']
# -

time = oope['time'].values + 1

pred = read_variable('predationTrend')
growth = read_variable('growthTrend')
zadv = read_variable('zadv_trend')
madv = read_variable('madv_trend')
zdiff = read_variable('zdiff_trend')
mdiff = read_variable('mdiff_trend')
# gamma1 = read_variable('gamma1')
# mort = read_variable('mort_day')
# u_act = read_variable('u_active')
# v_act = read_variable('v_active')
# u_pas = read_variable('u_passive')
# v_pas = read_variable('v_passive')
# repfonct = read_variable('repfonct_day')

def get_clim(var):
    var = var[~np.isnan(var)]
    perc = np.percentile(np.abs(var), 99)
    return -perc, perc


def plot(ax, toplot, wstep, contour=True, levels=None, clim=None, trend=True):
    toplot = toplot * wstep
    if clim is None:
        cmin, cmax = get_clim(toplot)
    else:
        cmin, cmax = [-clim, clim]
    if levels is None:
        levels = np.linspace(cmin, cmax, 11)
        
    levels = levels[levels != 0]

    time = np.arange(toplot.shape[0]) + 1

    cs = ax.pcolormesh(lon, time, toplot, shading='auto')
    #cs = ax.imshow(toplot, extent=[lon.min(), lon.max(), time.min(), time.max()], interpolation='none')
    cl = None
    cs.set_clim(cmin, cmax)
    if contour:
        cl = ax.contour(lon, time, toplot, levels=levels, colors='k', linewidths=2)
        #cl2 = ax.contour(lon, time, toplot, levels=0, linewidths=2, colors='k')
    stride = 3
    ax.xaxis.set_major_formatter(formatter0)
    ax.grid(True, linewidth=0.5, color='gray', linestyle='--')
    labels = ['180', '-150', '-120', '-90']
    xticks = np.array([float(l) for l in labels])
    xticks[xticks < 0] += 360
    ax.set_xticks(xticks)
    plt.setp(ax.get_xticklabels(), ha='right', rotation=45)
    ax.set_ylim(time.min(), time.max())
    ax.set_yticks(time[::3])
    ax.set_yticklabels(tlabels[::3], va='top', rotation=45)
    return cs, cl


dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
textprop = {}
textprop['bbox'] = dicttext
textprop['ha'] = 'center'
textprop['va'] = 'center'
textprop['zorder'] = 1000


def read_pisces_variable(varname, zmax=50, anom=False):

    e3t = mesh['e3t_0'] * mesh['tmask']
    e3t = e3t.where(e3t['olevel'] <= zmax)

    data = xr.open_dataset('../fig4/equatorial_full_%s.nc' %varname)[varname]
    data['olevel'] = mesh['olevel']
    data = (data * e3t).sum(dim='olevel') / (e3t.sum(dim='olevel'))
    if(anom):
        clim = data.sel(time_counter=slice('1971-01-01', '2000-12-31')).groupby('time_counter.month').mean(dim='time_counter')
        data = data.groupby('time_counter.month') - clim
    ylist = [1982, 1997, 2015]
    for year in ylist:
        datestart = '%s-01-01' %(year)
        dateend = '%s-12-31' %(year + 1)
        temp = data.sel(time_counter=slice(datestart, dateend))
        if year == 1982:
            thetao = temp.values
        else:
            thetao += temp.values
    thetao /= len(ylist)
    return thetao


months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
tlabels = ['%s-%d' %(m, y) for y in [0, 1] for m in months]
tlabels


def make_levels(ccc, step):
    ccc = abs(ccc)
    levels = np.arange(-ccc, ccc + step, step)
    levels = levels[levels != 0]
    return levels


# +
fig = plt.figure(figsize=(13, 14), facecolor='white')
plt.rcParams['font.size'] = 15

thetao = read_pisces_variable('thetao', 50, True)
anom = True
plk = read_pisces_variable('PLK', 50, anom)

time1 = 3
lon1 = 255

time2 = len(time) - 2
lon2 = 255

axgr = ImageGrid(fig, 111, nrows_ncols=(3, 3), axes_pad=(0.7, 0.3), cbar_pad=0.1, direction='row', aspect=False, cbar_mode="each", share_all=True)

ccc = 150
cpt = -1

#wstep_l0 = float(wstep.sel(l=l0, method='nearest'))
cprop = {}
cprop['colors'] = 'k'
cprop['linewidths']= 1

clist = [150, 15, 6]
cont = True
iii = 0

levels0 = make_levels(125, 25)
levels0 = make_levels(-150, 50)
levels1 = make_levels(-15, 5)
levels2 = make_levels(-5, 2)

levels = [levels0, levels1, levels2]

index = 0
for l in [3, 20, 90]:
    ccc = levels[iii].max()
    cpt += 1    
    ax = axgr[cpt]
    toplot = (wstep * (zadv + madv + zdiff + mdiff + pred + growth)).cumsum(dim='time').isel(l=index)
    cs, cl = plot(ax, toplot.values, 1, clim=ccc, levels=levels[iii])
    cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
    cb.add_lines(cl)
    cb.set_ticks(levels[iii])
    #cb.set_label('J/m2')
    ax.text(lon1, time1, letters[cpt], **textprop)
    ax.text(lon2, time2, '%scm' %l, **textprop)
    if(l == 3):
        ax.set_title('TOT (J/m2)')
    
    cpt += 1
    ax = axgr[cpt]
    toplot = (wstep * (pred + growth)).cumsum(dim='time').isel(l=index)
    skipcont = (l == 90)
    cs, cl = plot(ax, toplot.values, 1, clim=ccc, levels=levels[iii], contour=(skipcont == False))
    cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
    ax.text(lon1, time1, letters[cpt], **textprop)
    ax.text(lon2, time2, '%scm' %l, **textprop)
    if(l == 3):
        ax.set_title('P+G (J/m2)')
    if not skipcont:
        cb.add_lines(cl)
    cb.set_ticks(levels[iii])
    
    cpt += 1
    ax = axgr[cpt]
    toplot = (wstep * (zadv + madv + zdiff + mdiff)).cumsum(dim='time').isel(l=index)
    cs, cl = plot(ax, toplot.values, 1, clim=ccc, contour=True, levels=levels[iii])
    cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
    cb.add_lines(cl)
    cb.set_ticks(levels[iii])
    ax.text(lon1, time1, letters[cpt], **textprop)
    ax.text(lon2, time2, '%scm' %l, **textprop)
    if(l == 3):
        ax.set_title('A+D (J/m2)')
    step = 0.2
    iii += 1
    index += 1

plt.savefig('gr6.jpg', bbox_inches='tight')