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
from cartopy.mpl.ticker import LongitudeFormatter

formatter0 = LongitudeFormatter(dateline_direction_label=True)

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

boolean = (lon0 >= 150)
surface = surface.where(boolean, drop=True)
lon0 = lon0.where(boolean, drop=True)

# ## Loading the NEMO/Pisces anomalies

data = xr.open_mfdataset('data/pacific_0-50_*_anom.nc', compat='override').where(boolean, drop=True)
data = data.sel(time_counter=slice('1997-01-01', '1998-12-31'))
data

data['PLK'] = data['GOC'] + data['PHY2'] + data['ZOO'] + data['ZOO2']
data

# ## Loading Apecosm anomalies

isizes = [14, 45, 80]

dataape = xr.open_dataset('data/pacific_OOPE_anom.nc').where(boolean, drop=True)
dataape = dataape.sel(time=slice('1997-01-01', '1998-12-31'))
dataape

oope = dataape['OOPE']
oope

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc').isel(wpred=isizes)

length = const['length'].values * 100
length

wstep = const['weight_step'].values
wstep

oope = oope.isel(w=isizes)
oope

oope = oope.rename({'w': 'length'})
oope['length'] = length
oope

# ## Plotting the hovmoller

x = lon0.values
x

dates = data['time_counter'].values
dates

surfsum = surface.sum(dim='y')

y = np.arange(len(dates))
datestr = ['%.4d-%.2d' %(d.year, d.month) for d in dates]
datestr

# +
fig = plt.figure(facecolor='white', figsize=(14, 10))
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(2, 3),  # creates 2x2 grid of axes
                 axes_pad=[1.2, 0.6],  # pad between axes in inch.
                 cbar_mode='each', aspect=False, cbar_pad=0.1)
cbar_axes = grid.cbar_axes
stride = 3

quant = 0.99

lontext = 292
lattext = y[-3]
fs = 20
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

cpt = 0
for v in ['thetao', 'PLK', 'uo']:
    ax = grid[cpt]
    temp = (data[v] * surface).sum(dim='y') / surfsum
    cmax = float(abs(temp).quantile(quant))
    cs = ax.pcolormesh(x, y, temp, shading='auto')
    cl = ax.contour(x, y, temp, levels=np.linspace(-cmax, cmax, 11), colors='k', linewidths=0.5)
    cl0 = ax.contour(x, y, temp, levels=0, colors='k', linewidths=1)
    cs.set_clim(-cmax, cmax)
    cb = plt.colorbar(cs, cbar_axes[cpt])
    ax.set_yticks(y[::stride])
    ax.set_yticklabels(datestr[::stride], rotation=45, va='top')
    ax.grid(True)
    ax.set_xlabel('Longitude')
    ax.set_title(names[v])
    cb.set_label(units[v])
    ax.text(lontext, lattext, letters[cpt] + ")", ha='right', va='center', bbox=dicttext, fontsize=fs)
    ax.xaxis.set_major_formatter(formatter0)
    cpt += 1


nlength = len(length)
for l in range(nlength):
    ax = grid[cpt]
    temp = (oope.isel(length=l) * surface).sum(dim='y') / surfsum
    temp *= wstep[l]
    cmax = float(abs(temp).quantile(quant))
    cs = ax.pcolormesh(x, y, temp, shading='auto')
    cl = ax.contour(x, y, temp, levels=np.linspace(-cmax, cmax, 11), colors='k', linewidths=0.5)
    cl0 = ax.contour(x, y, temp, levels=0, colors='k', linewidths=1)
    cs.set_clim(-cmax, cmax)
    cb = plt.colorbar(cs, cbar_axes[cpt])
    ax.set_title('Biomass dens., L=%.f cm' %length[l])
    cb.set_label(units['oope'])
    ax.set_yticks(y[::stride])
    ax.set_yticklabels(datestr[::stride], rotation=45, va='top')
    ax.grid(True)
    ax.set_xlabel('Longitude')
    ax.text(lontext, lattext, letters[cpt] + ")", ha='right', va='center', bbox=dicttext, fontsize=fs)
    ax.xaxis.set_major_formatter(formatter0)
    labels = ['150', '180', '-150', '-120', '-90', '-60']
    xticks = np.array([float(l) for l in labels])
    xticks[xticks < 0] += 360
    ax.set_xticks(xticks)
    plt.setp(ax.get_xticklabels(), ha='right', rotation=45)

    cpt += 1

plt.savefig('plot_all_hovmoller_phys_oope.png', bbox_inches='tight')
# -


