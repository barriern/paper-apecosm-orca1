# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
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

plt.rcParams['text.usetex'] = False

latmax = 1
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

# ## Loading the NEMO/Pisces anomalies

data = xr.open_mfdataset('data/pacific_0-50_*_anom.nc', compat='override')
data

data['PLK'] = data['GOC'] + data['PHY2'] + data['ZOO'] + data['ZOO2']
data

# ## Loading Apecosm anomalies

isizes = [14, 45, 80]

dataape = xr.open_dataset('data/pacific_OOPE_anom.nc')
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
fig = plt.figure(figsize=(14, 10))
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(2, 3),  # creates 2x2 grid of axes
                 axes_pad=[0.8, 0.3],  # pad between axes in inch.
                 cbar_mode='each', aspect=False, cbar_pad=0.1)
cbar_axes = grid.cbar_axes
stride = 12

quant = 0.99

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
    cb.set_label(v)
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
    cb.set_label('OOPE, L=%.f cm' %length[l])
    ax.set_yticks(y[::stride])
    ax.set_yticklabels(datestr[::stride], rotation=45, va='top')
    ax.grid(True)
    ax.set_xlabel('Longitude')
    cpt += 1

plt.savefig('plot_all_hovmoller_phys_oope.pdf', bbox_inches='tight')
# -


