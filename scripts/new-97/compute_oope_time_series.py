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

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = False

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
length = const['length'].values * 100
const = const.rename({'wpred' : 'l'})
const['l'] = length
const

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0)
mesh

mask = mesh['tmaskutil']
mask.plot()

lon = mesh['glamt']
lat = mesh['gphit']

weights = mesh['e1t'] * mesh['e2t'] * mesh['tmaskutil']
weights.plot()

west = (abs(lat) < 5) & (lon >= 0)
weights.where(west).plot()

east = (abs(lat) < 5) & (lon <= 0)
weights.where(east).plot()

plt.figure(facecolor='white')
mask.plot(cmap='gray')
ilat, ilon = np.nonzero((east.values == True) & (mask.values== 1))
nlat, nlon = mask.shape
x = np.arange(nlon)
y = np.arange(nlat)
plt.plot(x[ilon], y[ilat], marker='.', color='steelblue')
ilat, ilon = np.nonzero((west.values == True) & (mask.values== 1))
plt.plot(x[ilon], y[ilat], marker='.', color='firebrick')
plt.savefig('mask_domain.png', bbox_inches='tight')

data = xr.open_dataset('data/pacific_nino97_OOPE.nc', decode_times=True)
oope = data['OOPE']
oope

dataclim = xr.open_dataset('data/pacific_clim_OOPE.nc')
oopeclim = dataclim['OOPE']
oopeclim

oope = oope.groupby('time.month') - oopeclim
oope

oope = oope.rename({'w': 'l'})
oope['l'] = length
oope = oope * const['weight_step']

oope_west = (oope * weights.where(west)).sum(dim=['x', 'y'])
oope_west.name = 'OOPE_west'
#oope_west.sel(l=3, method='nearest').plot()

oope_east = (oope * weights.where(east)).sum(dim=['x', 'y'])
oope_east.name = 'OOPE_east'
#oope_east.sel(l=3, method='nearest').plot()

dsout = xr.Dataset()
dsout['oope_east'] = oope_east
dsout['oope_west'] = oope_west
dsout.to_netcdf('oope_timeseries.nc')
