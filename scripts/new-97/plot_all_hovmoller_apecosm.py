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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import string

plt.rcParams['text.usetex'] = False

latmax = 5

letters = list(string.ascii_lowercase)

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
# -

# ## Loading mesh file

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0)
mesh

lat = mesh['gphit']
lat

mesh = mesh.where(abs(lat) <= latmax)
mesh

lon0 = mesh['glamt'].mean(dim='y')
lon0 = (lon0 + 360) % 360
lon0

surface = mesh['e1t'] * mesh['e2t'] * mesh['tmask']
surface

# ## Loading Apecosm anomalies

isizes = [14, 45, 80]

data = xr.open_mfdataset('data/pacific_[a-zA-Z]*_anom.nc').isel(w=isizes).sel(time=slice('1997-01', '1998-12'))
data

data = data.where(lon0 >= 150, drop=True)
mesh = mesh.where(lon0 >= 150, drop=True)
surface = surface.where(lon0 >= 150, drop=True)

lon0 = lon0.where(lon0 >= 150, drop=True)

data['adv'] = data['madv_trend'] + data['zadv_trend']
data

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc').isel(wpred=isizes)
const = const.rename({'wpred' : 'w'})

length = const['length'].values * 100
length

wstep = const['weight_step'].values
wstep

data['int-adv'] = data['adv'].cumsum(dim='time') * 24 * 60 * 60 * 30 * const['weight_step']

data['OOPE'] *= const['weight_step']
data

# ## Plotting the hovmoller

x = lon0.values
x

dates = data['time'].values
dates

surfsum = surface.sum(dim='y')

y = np.arange(len(dates))
datestr = ['%.4d-%.2d' %(d.year, d.month) for d in dates]
datestr

# +
varnames = ['repfonct_day', 'gamma1', 'mort_day', 'int-adv']
nvars = len(varnames)

sizes = [0, 1, 2]
nsizes = len(sizes)

fig = plt.figure(facecolor='white', figsize=(14, 16))
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(nvars, nsizes),  # creates 2x2 grid of axes
                 axes_pad=[1.1, 0.4],  # pad between axes in inch.
                 cbar_mode='each', aspect=False, cbar_pad=0.1)
cbar_axes = grid.cbar_axes
stride = 3

from matplotlib.ticker import FuncFormatter
fmt = lambda x, pos: '%.e' %(x)

quant = 0.98

lontext = 292
lattext = y[-3]
lattext2 = y[2]
fs = 20
plt.rcParams['font.size'] = 15
plt.rcParams['image.cmap'] = 'RdBu_r'

units = {}
units['int-adv'] = 'J/m2'
units['gamma1'] = ''
units['mort_day'] = ''
units['repfonct_day'] = ''

name = {}
name['int-adv'] = 'Adv. trend'
name['gamma1'] = 'Growth rate'
name['mort_day'] = 'Mort. rate'
name['repfonct_day'] = 'Func. response'

cpt = 0
for v in varnames:
    print('Plotting variable ', v)
    tempoope1 = (data['OOPE'] * surface).sum(dim='y') / surfsum
    temp1 = (data[v] * surface).sum(dim='y') / surfsum
    for s in sizes:
        ax = grid[cpt]
        temp = temp1.isel(w=s)
        tempoope = tempoope1.isel(w=s)
        cmaxoope = float(abs(tempoope).quantile(0.99))
        cmax = float(abs(temp).quantile(quant))
        if(v == 'int-adv'):
            print('processing int adv')
            cmax = cmaxoope
            temp = temp + tempoope.isel(time=0)
        cs = ax.pcolormesh(x, y, temp, shading='auto')
        cl = ax.contour(x, y, tempoope, levels=np.linspace(-cmaxoope, cmaxoope, 11), colors='k', linewidths=0.5)
        cl0 = ax.contour(x, y, tempoope, levels=0, colors='k', linewidths=1)
        cs.set_clim(-cmax, cmax)
        cb = plt.colorbar(cs, cbar_axes[cpt])
        ax.set_yticks(y[::stride])
        ax.set_yticklabels(datestr[::stride], rotation=45, va='top')
        ax.grid(True)
        ax.set_xlabel('Longitude')
        ax.set_title('%s' %(name[v]))
        cb.set_label('%s' %(units[v]))
        cb.ax.yaxis.set_offset_position('right')
        ax.text(lontext, lattext, letters[cpt] + ")", ha='right', va='center', bbox=dicttext, fontsize=fs)
        ax.text(lontext, lattext2, '%.f cm' %length[s], ha='right', va='center', bbox=dicttext, fontsize=fs)
        cpt += 1
plt.savefig('plot_all_hovmoller_apecosm.png', bbox_inches='tight')