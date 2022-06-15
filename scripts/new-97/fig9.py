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

mesh = xr.open_dataset('data/equatorial_mesh_mask.nc')
mesh = mesh.rename({'z': 'olevel'})
lon = mesh['x'].values
ilon = np.nonzero((lon >= 150) & (lon <= -90 + 360))[0]
mesh = mesh.isel(x=ilon)
lon = mesh['x'].values
depth = mesh['gdept_1d'].isel(x=0)
depth
mesh['olevel'] = depth
mesh

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100
const

wstep = const['weight_step']
wstep


# +
def read_variable(varname):
    
    data = xr.open_dataset('data/nino_equatorial_composites_%s.nc' %varname)[varname].isel(x=ilon)
    data['l'] = const['l']
    return data

oope = read_variable('OOPE')
oope
# -

time = oope['time'].values 

pred = read_variable('predationTrend')
growth = read_variable('growthTrend')
zadv = read_variable('zadv_trend')
madv = read_variable('madv_trend')
zdiff = read_variable('zdiff_trend')
mdiff = read_variable('mdiff_trend')
gamma1 = read_variable('gamma1')
mort = read_variable('mort_day')
u_act = read_variable('u_active')
v_act = read_variable('v_active')
u_pas = read_variable('u_passive')
v_pas = read_variable('v_passive')
repfonct = read_variable('repfonct_day')


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

    time = np.arange(toplot.shape[0]) + 1

    cs = ax.pcolormesh(lon, time, toplot, shading='auto')
    #cs = ax.imshow(toplot, extent=[lon.min(), lon.max(), time.min(), time.max()], interpolation='none')
    cs.set_clim(cmin, cmax)
    if contour:
        cl = ax.contour(lon, time, toplot, levels=levels, colors='k', linewidths=1)
        cl2 = ax.contour(lon, time, toplot, levels=0, linewidths=2, colors='k')
    stride = 3
    ax.xaxis.set_major_formatter(formatter0)
    ax.grid(True, linewidth=0.5, color='gray', linestyle='--')
    labels = ['180', '-150', '-120', '-90']
    xticks = np.array([float(l) for l in labels])
    xticks[xticks < 0] += 360
    ax.set_xticks(xticks)
    plt.setp(ax.get_xticklabels(), ha='right', rotation=45)
    ax.set_ylim(time.min(), time.max())
#     ax.set_ylabel('Month')
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    tlabels = ['%s-%d' %(m, y) for y in [0, 1] for m in months]
    tlabels
    ax.set_yticks(time[::3])
    ax.set_yticklabels(tlabels[::3], va='top', rotation=45)
    return cs, None


dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
textprop = {}
textprop['bbox'] = dicttext
textprop['ha'] = 'center'
textprop['fontsize'] = 20
textprop['va'] = 'center'
textprop['zorder'] = 1000


def read_pisces_variable(varname, zmax=50, anom=False):

    e3t = mesh['e3t_0'] * mesh['tmask']
    e3t = e3t.where(e3t['olevel'] <= zmax)

    data = xr.open_dataset('data/equatorial_full_%s.nc' %varname)[varname]
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


# +
fig = plt.figure(figsize=(18, 13), facecolor='white')
plt.rcParams['font.size'] = 15
plt.subplots_adjust(wspace=0.1, hspace=0.08)

thetao = read_pisces_variable('thetao', 50, True)
anom = True
phy2 = read_pisces_variable('PHY2', 50, anom)
zoo2 = read_pisces_variable('ZOO2', 50, anom)
zoo = read_pisces_variable('ZOO', 50, anom)
goc = read_pisces_variable('GOC', 50, anom)
plk = phy2 + zoo2 + zoo + goc

time1 = 3
lon1 = 255

time2 = len(time) - 2
lon2 = 255

ccc = 150
cpt = 0

cprop = {}
cprop['colors'] = 'k'
cprop['linewidths']= 1

clist = [150, 15, 5]
cont = True
iii = 0

l = 20

ccc = None

def set_ticks(ax, cpt):
    if cpt not in [1, 4]:
        plt.setp(ax.get_yticklabels(), visible=False)
        plt.ylabel('')
    else:
        plt.ylabel('Month')
        
    if cpt not in [4, 5, 6]:
        plt.setp(ax.get_xticklabels(), visible=False)
        

ccc = 3
cpt += 1 
ax = plt.subplot(2, 3, cpt)
toplot = (wstep * (growth + pred)).sel(l=l, method='nearest')
cs, cl = plot(ax, toplot.values, 1, clim=ccc, contour=0)
cb = plt.colorbar(cs, )
ax.text(lon1, time1, letters[cpt], **textprop)
ax.set_title('Pred. + Growth ($J.m^{-2}.s^{-1}$)')
cl = ax.contour(lon, time + 1, toplot.cumsum(dim='time'), colors='k', levels=np.linspace(-15, 15, 11), linewidths=1)
ax.contour(lon, time + 1, toplot.cumsum(dim='time'), colors='k', levels=[0], linewidths=2)
plt.clabel(cl)
set_ticks(ax, cpt)
        
ccc = 3
cpt += 1 
ax = plt.subplot(2, 3, cpt)
toplot = (wstep * (growth)).sel(l=l, method='nearest')
cs, cl = plot(ax, toplot.values, 1, clim=ccc, contour=0)
cb = plt.colorbar(cs, )
ax.text(lon1, time1, letters[cpt], **textprop)
ax.set_title('Growth ($J.m^{-2}.s^{-1}$)')
cl = ax.contour(lon, time + 1, toplot.cumsum(dim='time'), colors='k', levels=np.linspace(-15, 15, 11), linewidths=1)
ax.contour(lon, time + 1, toplot.cumsum(dim='time'), colors='k', levels=[0], linewidths=2)
plt.clabel(cl)
set_ticks(ax, cpt)

cpt += 1
ax = plt.subplot(2, 3, cpt)
toplot = (wstep * (pred)).sel(l=l, method='nearest')
cs, cl = plot(ax, toplot.values, 1, clim=ccc, contour=False)
cb = plt.colorbar(cs, )
ax.text(lon1, time1, letters[cpt], **textprop)
ax.set_title('Pred ($J.m^{-2}.s^{-1}$)')
cl = ax.contour(lon, time + 1, toplot.cumsum(dim='time'), colors='k', levels=np.linspace(-15, 15, 11), linewidths=1)
ax.contour(lon, time + 1, toplot.cumsum(dim='time'), colors='k', levels=[0], linewidths=2)
plt.clabel(cl)
set_ticks(ax, cpt)

cpt += 1
ax = plt.subplot(2, 3, cpt)
toplot = (repfonct).sel(l=l, method='nearest')
cs, cl = plot(ax, toplot.values, 1, clim=None, contour=False)
cb = plt.colorbar(cs, )
ax.text(lon1, time1, letters[cpt], **textprop)
ax.set_title('Func. response')
toplot = (oope * wstep).sel(l=3, method='nearest')
cl = ax.contour(lon, time + 1, toplot, 6, colors='k', linewidths=1)
ax.contour(lon, time + 1, toplot, colors='k', levels=[0], linewidths=2)
plt.clabel(cl)
set_ticks(ax, cpt)

cpt += 1
ax = plt.subplot(2, 3, cpt)
toplot = (gamma1).sel(l=l, method='nearest')
cs, cl = plot(ax, toplot.values, 1, clim=None, contour=False)
cb = plt.colorbar(cs, )
ax.text(lon1, time1, letters[cpt], **textprop)
ax.set_title('Growth rate ($kg.day^{-1}$)')
toplot = thetao
cl = ax.contour(lon, time + 1, toplot, 6, colors='k', linewidths=1)
ax.contour(lon, time + 1, toplot, colors='k', levels=[0], linewidths=2)
plt.clabel(cl)
set_ticks(ax, cpt)

cpt += 1
ax = plt.subplot(2, 3, cpt)
toplot = (mort).sel(l=l, method='nearest')
cs, cl = plot(ax, toplot.values, 1, clim=None, contour=False)
cb = plt.colorbar(cs, )
ax.text(lon1, time1, letters[cpt], **textprop)
ax.set_title('Pred. mort. rate ($day^{-1}$)')
toplot = (oope * wstep).sel(l=90, method='nearest')
cl = ax.contour(lon, time + 1, toplot, 6, colors='k', linewidths=1)
ax.contour(lon, time + 1, toplot, colors='k', levels=[0], linewidths=2)
plt.clabel(cl)
set_ticks(ax, cpt)

plt.savefig('fig9.png', bbox_inches='tight')
# -

