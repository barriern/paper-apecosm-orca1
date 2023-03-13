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
from dask.diagnostics import ProgressBar, visualize, Profiler, ResourceProfiler, CacheProfiler

varname = 'OOPE'
# -

const = xr.open_dataset('../data/static/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100
const

filename = os.path.join('../data/apecosm/*OOPE_Y*')
filename

data = xr.open_mfdataset(filename)
data = data.rename({'w': 'l'})
data['l'] = const['l']
data

var = data[varname]
var

var = var.sel(l=[3, 20, 90], method='nearest')
var

clim = var.groupby('time.month').mean(dim='time')
clim

fileout = 'clim_epi_%s.nc' %varname
fileout

delayed = clim.to_netcdf(fileout, compute=False)

with ProgressBar() as cbar:
    delayed.compute()


