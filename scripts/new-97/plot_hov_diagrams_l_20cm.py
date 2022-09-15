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

l0 = 20
# -

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
    ax.set_ylabel('Month')
    return cs, None


dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
textprop = {}
textprop['bbox'] = dicttext
textprop['ha'] = 'center'
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
fig = plt.figure(figsize=(13, 14), facecolor='white')
plt.rcParams['font.size'] = 15

thetao = read_pisces_variable('thetao', 50, True)
anom = True
phy2 = read_pisces_variable('PHY2', 50, anom)
zoo2 = read_pisces_variable('ZOO2', 50, anom)
zoo = read_pisces_variable('ZOO', 50, anom)
goc = read_pisces_variable('GOC', 50, anom)
plk = phy2 + zoo2 + zoo + goc

time1 = 3
lon1 = 255

time2 = len(time) - 3
lon2 = 255

axgr = ImageGrid(fig, 111, nrows_ncols=(2, 2), axes_pad=(1.1, 0.5), cbar_pad=0.1, direction='row', aspect=False, cbar_mode="each", share_all=True)

ccc = 15
cpt = -1

wstep_l0 = float(wstep.sel(l=l0, method='nearest'))
cprop = {}
cprop['colors'] = 'k'
cprop['linewidths']= 1

cont = True


cpt += 1
toplot = (growth).sel(l=l0, method='nearest').cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep_l0, clim=ccc)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'G', **textprop)
cb.set_label('J/m2')
ax.text(lon1, time1, letters[cpt], **textprop)

cpt += 1
toplot = (pred).sel(l=l0, method='nearest').cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep_l0, clim=ccc)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'P', **textprop)
cb.set_label('J/m2')
ax.text(lon1, time1, letters[cpt], **textprop)

cpt += 1
toplot = (pred + growth).sel(l=l0, method='nearest').cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep_l0, clim=ccc)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'P+G', **textprop)
cb.set_label('J/m2')
ax.text(lon1, time1, letters[cpt], **textprop)

cpt += 1
toplot = (zadv + madv + zdiff + mdiff).sel(l=l0, method='nearest').cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep_l0, clim=ccc, contour=cont)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'A+D', **textprop)
ax.text(lon1, time1, letters[cpt], **textprop)
cb.set_label('J/m2')

plt.savefig('hov_compo_l_%d.png' %l0, bbox_inches='tight')
# +
fig = plt.figure(figsize=(13, 14), facecolor='white')
plt.rcParams['font.size'] = 15

thetao = read_pisces_variable('thetao', 50, True)
anom = False
phy2 = read_pisces_variable('PHY2', 50, anom)
zoo2 = read_pisces_variable('ZOO2', 50, anom)
zoo = read_pisces_variable('ZOO', 50, anom)
goc = read_pisces_variable('GOC', 50, anom)
plk = phy2 + zoo2 + zoo + goc

cont = True

time1 = 3
lon1 = 255

time2 = len(time) - 3
lon2 = 255

cont = False

axgr = ImageGrid(fig, 111, nrows_ncols=(3, 2), axes_pad=(1.1, 0.5), cbar_pad=0.1, direction='row', aspect=False, cbar_mode="each", share_all=True)

ccc = 15
cpt = -1

wstep_l0 = float(wstep.sel(l=l0, method='nearest'))
cprop = {}
cprop['colors'] = 'k'
cprop['linewidths']= 1

ccc = -0.6
ccc = 15

cpt += 1
toplot = (growth).sel(l=l0, method='nearest').cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep_l0, clim=ccc)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'G', **textprop)
cb.set_label('J/m2')
ax.text(lon1, time1, letters[cpt], **textprop)

cpt += 1
toplot = (pred).sel(l=l0, method='nearest').cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep_l0, clim=ccc)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'P', **textprop)
cb.set_label('J/m2')
ax.text(lon1, time1, letters[cpt], **textprop)

cpt += 1
toplot = (pred + growth).sel(l=l0, method='nearest').cumsum(dim='time').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, wstep_l0, clim=ccc)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'P+G', **textprop)
cb.set_label('J/m2')
ax.text(lon1, time1, letters[cpt], **textprop)

ccc = None
cpt += 1
toplot = (repfonct).sel(l=l0, method='nearest').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, 1, clim=ccc, contour=cont)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'repfonct', **textprop)
cb.set_label('')
ax.text(lon1, time1, letters[cpt], **textprop)
toplot = (oope * wstep).sel(l=3, method='nearest').values
cl = ax.contour(lon[:], time + 1, toplot[:, :], 6, **cprop)
cl2 = ax.contour(lon[:], time + 1, toplot[:, :], levels=[0], linewidths=2, colors=cprop['colors'])
plt.clabel(cl)

cpt += 1
toplot = (gamma1).sel(l=l0, method='nearest').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, 1, clim=ccc, contour=cont)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'Gamma', **textprop)
cb.set_label('')
ax.text(lon1, time1, letters[cpt], **textprop)
toplot = thetao
cl = ax.contour(lon[:], time + 1, toplot[:, :], 6, **cprop)
cl2 = ax.contour(lon[:], time + 1, toplot[:, :], levels=[0], linewidths=2, colors=cprop['colors'])
plt.clabel(cl)

cpt += 1
toplot = (mort).sel(l=l0, method='nearest').values
ax = axgr[cpt]
cs, cl = plot(ax, toplot, 1, clim=ccc, contour=cont)
cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
ax.text(lon2, time2, 'Mort', **textprop)
cb.set_label('')
ax.text(lon1, time1, letters[cpt], **textprop)
toplot = (oope * wstep).sel(l=90, method='nearest').values
cl = ax.contour(lon[:], time + 1, toplot[:, :], 6, **cprop)
cl2 = ax.contour(lon[:], time + 1, toplot[:, :], levels=[0], linewidths=2, colors=cprop['colors'])
plt.clabel(cl)

# cpt += 1
# toplot = (v_pas).sel(l=l0, method='nearest').values
# ax = axgr[cpt]
# cs, cl = plot(ax, toplot, 1, clim=ccc, contour=cont)
# cb = plt.colorbar(cs, cax=axgr.cbar_axes[cpt])
# ax.text(lon2, time2, 'V pas', **textprop)
# cb.set_label('m/s')
# ax.text(lon1, time1, letters[cpt], **textprop)

# -

vpas_max = abs(v_pas).max(dim=['time', 'x'])
vact_max = abs(v_act).max(dim=['time', 'x'])
(vpas_max / vact_max).sel(l=l0, method='nearest')

upas_max = abs(u_pas).max(dim=['time', 'x'])
uact_max = abs(u_act).max(dim=['time', 'x'])
(upas_max / uact_max).sel(l=l0, method='nearest')


