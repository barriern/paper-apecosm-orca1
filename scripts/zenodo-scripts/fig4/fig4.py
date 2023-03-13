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
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import string
from cartopy.mpl.ticker import LongitudeFormatter

zmax = 50
formatter0 = LongitudeFormatter(dateline_direction_label=True)

latmax = 5

letters = list(string.ascii_lowercase)

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
# -

const = xr.open_dataset('../data/static/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100

# ## Loading mesh file

mesh = xr.open_dataset('../data/static/pacific_mesh_mask.nc').isel(y=61)
mesh = mesh.rename({'z': 'olevel'})
mesh['olevel'] = mesh['gdept_1d'].isel(x=0)
mesh = mesh.where(mesh['olevel'] <= zmax)

e3t = mesh['e3t_0']
e3t

lon0 = mesh['glamt'].isel(olevel=0)
lon0 = (lon0 + 360) % 360
lon0 = lon0.drop(['olevel'])
lon0

boolean = ((lon0 >= 150) & (lon0 <= 270))
lon0 = lon0.where(boolean, drop=True)
e3t = e3t.where(boolean, drop=True)

# ## Loading the NEMO/Pisces anomalies

data = xr.open_mfdataset('nino_equatorial_composites_*.nc', compat='override')
data

data = data.where(boolean, drop=True)
data

oope = data['OOPE']
oope['l'] = const['l']
oope = oope * const['weight_step']

lengths = [3, 20, 90]
oope = oope.sel(l=lengths, method='nearest')
oope

e3t['olevel'] = data['olevel']

# ## Computing warm pool displacement

# First, computing the climatology.

varname = 'thetao'
datatemp = xr.open_dataset('equatorial_full_%s.nc' %varname)[varname].where(boolean==True, drop=True)
temp = (datatemp * e3t).sum(dim='olevel') / e3t.sum(dim='olevel')
tempclim = temp.sel(time_counter=slice('1971-01-01', '2000-12-31')).groupby('time_counter.month').mean(dim='time_counter')
tempclim

tempanom = (data['thetao'] * e3t).sum(dim='olevel') / e3t.sum(dim='olevel')
warm_pool = tempanom.values
warm_pool[:12] = warm_pool[:12] + tempclim.values
warm_pool[12:] = warm_pool[12:] + tempclim.values

# ## Plotting the hovmoller

x = lon0.values
x

months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
tlabels = ['%s-%d' %(m, y) for y in [0, 1] for m in months]
tlabels


def make_levels(ccc, step):
    
    levels = np.arange(-ccc, ccc + step, step)
    levels = levels[levels != 0]
    return levels
# +
fig = plt.figure(facecolor='white', figsize=(14, 10))
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(2, 3),  # creates 2x2 grid of axes
                 axes_pad=[1.2, 0.6],  # pad between axes in inch.
                 cbar_mode='each', aspect=False, cbar_pad=0.1)
cbar_axes = grid.cbar_axes
stride = 3
plt.rcParams['text.usetex'] = False
plt.rcParams['font.size'] = 15

y = np.arange(1, 25)

quant = 0.99

lontext = 260
lattext = y[-3]
fs = 15
plt.rcParams['font.size'] = 15
plt.rcParams['image.cmap'] = 'RdBu_r'

units = {}
units['thetao'] = 'C'
units['PLK'] = 'mmol/m3'
units['uo'] = 'm/s'
units['oope'] = 'J/m2'



names = {}
names['thetao'] = 'Temperature'
names['PLK'] = 'Plankton conc.'
names['uo'] = 'Zonal vel.'

ccc = {}
ccc['thetao'] = make_levels(6, 2)
ccc['uo'] = make_levels(1, 0.25)
ccc['PLK'] = make_levels(1.5, 0.5)

ccc_oope = []
ccc_oope.append(make_levels(150, 50))
ccc_oope.append(make_levels(15, 5))
ccc_oope.append(make_levels(5, 2))


cpt = 0
for v in ['thetao', 'uo', 'PLK']:
    ax = grid[cpt]
    temp = (data[v] * e3t).sum(dim='olevel') / e3t.sum(dim='olevel')
    if(v == 'thetao'):
        print('lala')
        ax.contour(x, y, warm_pool, levels=[28], colors='r', zorder=1000)
    levels = ccc[v]
    cmax = float(abs(temp).quantile(quant))
    cs = ax.pcolormesh(x, y, temp, shading='auto')
    cl = ax.contour(x, y, temp, levels=levels, colors='k', linewidths=2)
    #cl0 = ax.contour(x, y, temp, levels=0, colors='k', linewidths=1)
    cs.set_clim(levels.min(), levels.max())
    cb = plt.colorbar(cs, cbar_axes[cpt])
    cb.add_lines(cl)
    cb.set_ticks(levels)
    ax.grid(True)
    ax.set_xlabel('Longitude')
    #ax.set_ylabel('Months')
    ax.set_title(names[v])
    cb.set_label(units[v])
    ax.text(lontext, lattext, letters[cpt] + ")", ha='right', va='center', bbox=dicttext, fontsize=fs)
    ax.xaxis.set_major_formatter(formatter0)
    ax.set_ylim(y.min(), y.max())
    ax.set_yticks(y[::3])
    ax.set_yticklabels(tlabels[::3], va='top', rotation=45)
    cpt += 1


nlength = len(lengths)
for l in range(nlength):
    ax = grid[cpt]
    temp = oope.isel(l=l)
    cmax = float(abs(temp).quantile(quant))
    cs = ax.pcolormesh(x, y, temp, shading='auto')
    levels = ccc_oope[l]
    cl = ax.contour(x, y, temp, levels=levels, colors='k', linewidths=2)
    #cl0 = ax.contour(x, y, temp, levels=0, colors='k', linewidths=1)
    cs.set_clim(levels.min(), levels.max())
    cb = plt.colorbar(cs, cbar_axes[cpt])
    ax.set_title('Biomass dens., L=%.f cm' %lengths[l])
    cb.set_label(units['oope'])
    cb.add_lines(cl)
    cb.set_ticks(levels)
    ax.grid(True)
    ax.set_xlabel('Longitude')
    #ax.set_ylabel('Months')
    ax.text(lontext, lattext, letters[cpt] + ")", ha='right', va='center', bbox=dicttext, fontsize=fs)
    ax.xaxis.set_major_formatter(formatter0)
    labels = ['180', '-150', '-120', '-90']
    xticks = np.array([float(l) for l in labels])
    xticks[xticks < 0] += 360
    ax.set_xticks(xticks)
    plt.setp(ax.get_xticklabels(), ha='right', rotation=45)
    ax.set_ylim(y.min(), y.max())
    ax.set_yticks(y[::3])
    ax.set_yticklabels(tlabels[::3], va='top', rotation=45)
    cpt += 1

plt.savefig('gr4.jpg', bbox_inches='tight')