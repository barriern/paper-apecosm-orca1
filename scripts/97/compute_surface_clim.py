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

# # Compute surface anomalies
#
# ## Imports

# +
import xarray as xr
import os
from glob import glob
from dask.diagnostics import visualize, ProgressBar, Profiler, ResourceProfiler, CacheProfiler

varname = 'O2'
grid = 'ptrc_T'
# -

filelist = glob('data/*%s*nc' %grid)
filelist.sort()
filelist

data = xr.open_mfdataset(filelist).isel(olevel=0, )
data

var = data[varname]
var

clim = var.groupby('time_counter.month').mean(dim='time_counter')
clim

delayed = clim.to_netcdf('clim_%s.nc' %varname, compute=False)

with ProgressBar(),  Profiler() as prof, ResourceProfiler(dt=0.25) as rprof, CacheProfiler() as cprof:
    delayed.compute()

visualize([prof, rprof, cprof])

clim.data.visualize()


