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
#     display_name: Python 3 (ipykernel)
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

zmax = 50
formatter0 = LongitudeFormatter(dateline_direction_label=True)

plt.rcParams['text.usetex'] = False

latmax = 5

letters = list(string.ascii_lowercase)

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
# -

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100

# ## Loading mesh file

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(y=61)
mesh = mesh.rename({'z': 'olevel'})
mesh['olevel'] = mesh['gdept_1d'].isel(x=0)
mesh = mesh.where(mesh['olevel'] <= zmax)

e3t = mesh['e3t_0']

lon0 = mesh['glamt'].isel(olevel=0)
lon0 = (lon0 + 360) % 360
lon0

boolean = ((lon0 >= 150) & (lon0 <= 270))
lon0 = lon0.where(boolean, drop=True)
e3t = e3t.where(boolean, drop=True)

# ## Loading the NEMO/Pisces anomalies

data = xr.open_mfdataset('data/nino_equatorial_composites_*.nc', compat='override').where(boolean, drop=True)
data

data['PLK'] = data['GOC'] + data['PHY2'] + data['ZOO'] + data['ZOO2']
data

oope = data['OOPE']
oope['l'] = const['l']
oope = oope * const['weight_step']

lengths = [3, 20, 90]
oope = oope.sel(l=lengths, method='nearest')
oope

e3t['olevel'] = data['olevel']

# ## Plotting the hovmoller

x = lon0.values
x

# +
fig = plt.figure(facecolor='white', figsize=(14, 10))
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(2, 3),  # creates 2x2 grid of axes
                 axes_pad=[1.2, 0.6],  # pad between axes in inch.
                 cbar_mode='each', aspect=False, cbar_pad=0.1)
cbar_axes = grid.cbar_axes
stride = 3
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

cpt = 0
for v in ['thetao', 'PLK', 'uo']:
    ax = grid[cpt]
    temp = (data[v] * e3t).sum(dim='olevel') / e3t.sum(dim='olevel')
    cmax = float(abs(temp).quantile(quant))
    cs = ax.pcolormesh(x, y, temp, shading='auto')
    cl = ax.contour(x, y, temp, levels=np.linspace(-cmax, cmax, 11), colors='k', linewidths=0.5)
    cl0 = ax.contour(x, y, temp, levels=0, colors='k', linewidths=1)
    cs.set_clim(-cmax, cmax)
    cb = plt.colorbar(cs, cbar_axes[cpt])
    ax.grid(True)
    ax.set_xlabel('Longitude')
    ax.set_title(names[v])
    cb.set_label(units[v])
    ax.text(lontext, lattext, letters[cpt] + ")", ha='right', va='center', bbox=dicttext, fontsize=fs)
    ax.xaxis.set_major_formatter(formatter0)
    ax.set_ylim(y.min(), y.max())
    cpt += 1


nlength = len(lengths)
for l in range(nlength):
    ax = grid[cpt]
    temp = oope.isel(l=l)
    cmax = float(abs(temp).quantile(quant))
    cs = ax.pcolormesh(x, y, temp, shading='auto')
    cl = ax.contour(x, y, temp, levels=np.linspace(-cmax, cmax, 11), colors='k', linewidths=0.5)
    cl0 = ax.contour(x, y, temp, levels=0, colors='k', linewidths=1)
    cs.set_clim(-cmax, cmax)
    cb = plt.colorbar(cs, cbar_axes[cpt])
    ax.set_title('Biomass dens., L=%.f cm' %lengths[l])
    cb.set_label(units['oope'])
    ax.grid(True)
    ax.set_xlabel('Longitude')
    ax.text(lontext, lattext, letters[cpt] + ")", ha='right', va='center', bbox=dicttext, fontsize=fs)
    ax.xaxis.set_major_formatter(formatter0)
    labels = ['180', '-150', '-120', '-90']
    xticks = np.array([float(l) for l in labels])
    xticks[xticks < 0] += 360
    ax.set_xticks(xticks)
    plt.setp(ax.get_xticklabels(), ha='right', rotation=45)
    ax.set_ylim(y.min(), y.max())
    cpt += 1

plt.savefig('plot_all_hovmoller_phys_oope.png', bbox_inches='tight')
# -


