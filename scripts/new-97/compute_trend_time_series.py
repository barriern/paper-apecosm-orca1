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

import xarray as xr
import matplotlib.pyplot as plt
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

# Now processing the trends

varnames = ['growthTrend',
 'madv_trend',
 'mdiff_trend',
 'predationTrend',
 'zadv_trend',
 'zdiff_trend']
varnames

for v in varnames:
    
    print('Processing ', v)
    
    data = xr.open_mfdataset('data/*nino97*%s.nc' %v, decode_times=True)
    data = data[v]
    data

    dataclim = xr.open_mfdataset('data/*clim*%s.nc' %v, decode_times=True)
    dataclim = dataclim[v]
    dataclim

    data = (data.groupby('time.month') - dataclim).cumsum(dim='time')
    data

    data = data.rename({'w': 'l'})
    data['l'] = length
    data

    dsout = xr.Dataset()

    data_west = (data * weights.where(west)).sum(dim=['x', 'y'])
    data_west.name = '%s_west' %v
    dsout['%s_west' %v] = data_west

    data_east = (data * weights.where(east)).sum(dim=['x', 'y'])
    data_east.name = '%s_east' %v
    dsout['%s_east' %v] = data_east

    dsout.to_netcdf('%s_timeseries.nc' %v)
