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
#     display_name: Python [conda env:nbarrier] *
#     language: python
#     name: conda-env-nbarrier-py
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import string
from dask.diagnostics import ProgressBar
plt.rcParams['text.usetex'] = False

itime = slice('1997-01-01', '1998-12-31')

letters = list(string.ascii_lowercase)

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
# -

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l', 'gpred': 'c'})
const['l'] = const['length'].values * 100
const

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0)
mesh

lat = mesh['gphit'].to_masked_array()
np.abs(lat).min()
ilat, iok = np.nonzero(np.abs(lat) == 0)
ilat = ilat[0]
ilat

mesh = mesh.isel(y=ilat)
mesh

data = xr.open_mfdataset('data/pacific_OOPE_anom.nc').isel(y=ilat).sel(time=itime)
data = data.rename({'w': 'l'})
data['l'] = const['length'].values * 100
data = data['OOPE']
data['x'] = mesh['glamt'].values
data.name = 'OOPE'
data.to_netcdf('data/hovmoller_plot_data_oope.nc')
data

varnames = [
    'growthTrend',
    'madv_trend',
    'mdiff_trend',
    'predationTrend',
    'zadv_trend',
    'zdiff_trend'
]
varnames

for v in varnames:

    print('Processing variable ', v)
    data = xr.open_dataset('data/pacific_nino97_%s.nc' %v).isel(y=ilat).sel(time=itime)
    temp = data[v]
#     temp = temp.cumsum(dim='time')
#     temp['x'] = mesh['glamt'].values
#     temp = temp.rename({'w': 'l'})
#     temp['l'] = const['length'].values * 100
#     temp.to_netcdf('data/hovmoller_plot_data_%s.nc' %v)

    dataclim = xr.open_dataset('data/pacific_clim_%s.nc' %v).isel(y=ilat)
    tempclim = dataclim[v]
    tempclim

    tempanom = temp.groupby('time.month') - tempclim
    tempanom = tempanom.rename({'w' : 'l'})
    tempanom['l'] = const['length'].values * 100
    tempanom

    newtrend = tempanom.cumsum(dim='time')
    newtrend.name = v
    newtrend['x'] = mesh['glamt'].values

    filename = 'data/hovmoller_plot_data_%s.nc' %v
    with ProgressBar():
        newtrend.to_netcdf(filename)






