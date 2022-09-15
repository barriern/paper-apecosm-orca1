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


def get_clim(var):
    var = var[~np.isnan(var)]
    perc = np.percentile(np.abs(var), 99)
    return -perc, perc


# +
def read_variable(varname):
    
    data = xr.open_dataset('data/nino_equatorial_composites_%s.nc' %varname)[varname].isel(x=ilon)
    data['l'] = const['l']
    return data

oope = read_variable('OOPE')
oope
# -

oope_clim = xr.open_dataset('data/pacific_clim_OOPE.nc')['OOPE'].isel(x=ilon, y=61)
oope_clim = oope_clim.rename({'w': 'l'})
oope_clim['l'] = const['l']
oope_clim

time = oope['time'].values 

gamma1 = read_variable('gamma1')
x = gamma1['x'].values
length = gamma1['l'].values
gamma1

l = 20
dist = (length - l)**2
iok = np.argmin(dist)
length[iok]

# +
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
tp1 = gamma1.values[:, :, iok]
cs = plt.pcolormesh(x, time, tp1)
cb = plt.colorbar(cs)
plt.title('gamma1')
cs.set_clim(get_clim(tp1))

plt.subplot(2, 2, 2)
tp2 = tp1.copy()
tp2[:12] = oope_clim.values[..., iok]
tp2[12:] = tp2[:12]

cs = plt.pcolormesh(x, time, tp2)
cb = plt.colorbar(cs)
plt.title('oope')
cs.set_clim(get_clim(tp2))

plt.subplot(2, 2, 3)
tp3 = tp1 * tp2
cs = plt.pcolormesh(x, time, tp3)
cb = plt.colorbar(cs)
plt.title('gamma1 x oope')
cs.set_clim(get_clim(tp3))
# -


