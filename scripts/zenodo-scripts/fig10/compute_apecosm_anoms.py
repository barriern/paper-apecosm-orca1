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
#     display_name: Python [conda env:nbarrier]
#     language: python
#     name: conda-env-nbarrier-py
# ---

# +
import xarray as xr
import os
import numpy as np
from dask.diagnostics import ProgressBar

varname = 'OOPE'
# -

const = xr.open_dataset('../data/static/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100
const

filename = os.path.join('../data/apecosm/*OOPE_Y*')
filename

climname = 'clim_epi_OOPE.nc'
climname

data = xr.open_mfdataset(filename, chunks={'time': 24, 'w': 100, 'x': 100, 'y': 100})
data = data.rename({'w': 'l'})
data['l'] = const['l']
data = data.sel(l=[3, 20, 90], method='nearest')
data = data['OOPE']
data

clim = xr.open_dataset(climname, chunks={'month': -1, 'x': 100, 'y': 100})
clim[varname]

anom = data.groupby('time.month') - clim
anom[varname]

fileout = os.path.join('/home1/scratch/nbarrier/', 'anom_epi_%s.nc' %varname)
fileout

delayed = anom.to_netcdf(fileout, compute=False)

with ProgressBar():
    delayed.compute()


