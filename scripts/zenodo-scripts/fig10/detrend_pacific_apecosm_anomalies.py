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
from glob import glob
import numpy as np
import scipy.signal as sig
from dask.diagnostics import ProgressBar

varname = 'OOPE'
dirin = os.getenv('SCRATCH')
dirin
# -

filelist = glob(os.path.join(dirin, 'anom*%s*.nc' %varname))
filelist.sort()
filelist

data = xr.open_mfdataset(filelist, combine='by_coords')
data

var = data[varname]
var

var = var.chunk({'time' : -1, 'y': 20, 'x': 20})
var


# ## Creating all that is needed to detrending in par.

def gufunc_detrend(x):
    x[np.isnan(x)] = 0
    return sig.detrend(x)


def xarray_detrend(x, dim):
    return xr.apply_ufunc(
        gufunc_detrend,
        x,
        input_core_dims=[[dim]],
        output_core_dims=[[dim]],
        dask="parallelized",
        output_dtypes=[np.float32],
    )



detrended = xarray_detrend(var, dim='time')
detrended

detrended

# ## Saving output file

fileout = '%s/detrended_pacific_%s_anomalies.nc' %(dirin, varname)
fileout

with ProgressBar():
    detrended.T.to_netcdf(fileout)


