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
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar

itime = slice('1997-10-01', '1997-12-31')
chunk = {'x' : 50, 'y' : 50}
# -

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0)
lonf = mesh['glamf'].values
latf = mesh['gphif'].values

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l', 'gpred': 'c'})
length = const['length'] * 100
const['l'] = length
weight_step = const['weight_step']
weight_step

data = xr.open_dataset('data/pacific_nino97_OOPE.nc').chunk(chunk)
oope = data['OOPE']
oope

dataclim = xr.open_dataset('data/pacific_clim_OOPE.nc').chunk(chunk)
oopeclim = dataclim['OOPE']
oopeclim

oopeanom = oope.groupby('time.month') - oopeclim
oopeanom = oopeanom.rename({'w' : 'l'})
oopeanom['l'] = length
oopeanom

varnames = [
    'growthTrend',
    'madv_trend',
    'mdiff_trend',
    'predationTrend',
    'zadv_trend',
    'zdiff_trend'
]
varnames

for v in varnames[:]:
    
    print('Processing ', v)
    
    data = xr.open_dataset('data/pacific_nino97_%s.nc' %v).chunk(chunk)
    temp = data[v]
    temp

    dataclim = xr.open_dataset('data/pacific_clim_%s.nc' %v).chunk(chunk)
    tempclim = dataclim[v]
    tempclim
    
    tempanom = temp.groupby('time.month') - tempclim
    tempanom = tempanom.rename({'w' : 'l'})
    tempanom['l'] = length

    newtrend = tempanom.cumsum(dim='time')
    
    filename = 'data/integrated_maps_%s.nc' %v
    with ProgressBar():
        newtrend.to_netcdf(filename)
