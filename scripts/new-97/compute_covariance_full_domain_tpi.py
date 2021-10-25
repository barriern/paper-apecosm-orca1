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
import os
import scipy.signal as sig 
import numpy as np
from dask.diagnostics import ProgressBar
from glob import glob
import sys
import matplotlib.pyplot as plt

dirin = os.getenv('SCRATCH')

latmax = 10
lonmin = 150
lonmax = -100
# -

# ## Loading TPI index

tpi = xr.open_dataset('../data/filt_tpi.nc')
tpi = tpi['tpi_filt']

# ## Loading the datasets

dirin = os.getenv('SCRATCH')
data = xr.open_dataset(os.path.join(dirin, 'detrended_pacific_OOPE_anomalies.nc'))
data = data['OOPE']
data  # y, x, w, time

ny, nx, nw, ntime = data.shape

# Now, we correct the date of the TPI variable, to make sure everything works.

tpi['time'] = data['time']
tpi

# ## Preparing stuff for computation of lead-lag covariances

lags = sig.correlation_lags(ntime, ntime)
lags

ilags = np.nonzero(np.abs(lags) <= 10*12)[0]
nlags = len(ilags)
lags = lags[ilags]


def gufunc_cov(x, y):
    print('X.shape ', x.shape)  # (6, 40, 100, 1, 732)  => (x, y, w, 1, ntime)
    print('Y.shape ', y.shape)  #(100, 2, 732) => (weight, eof, time)
    nx = x.shape[0]
    ny = x.shape[1]
    nw = x.shape[2]
    ### Y = PCS
    ntime = x.shape[-1]
    lags = sig.correlation_lags(ntime, ntime)[ilags]
    nlags = len(lags)
    output = np.zeros((nx, ny, nw, nlags))
    for s in range(nx):
        for i in range(ny):
            for k in range(nw):
                xtemp = x[s, i, k, :]
                output[s, i, k, :] = sig.correlate(xtemp, y)[ilags] / ntime
    return output


def xarray_cov(x, y, dim):
    return xr.apply_ufunc(
        gufunc_cov,
        x,
        y,
        input_core_dims=[[dim], [dim]],
        output_core_dims=[['lags']],
        dask="parallelized",
        output_dtypes=[np.float32],
        dask_gufunc_kwargs = {'output_sizes' : {'lags': nlags}}
    )


data = data.chunk({'w' : 20, 'x' : 30, 'y' : 40})
data

cov = xarray_cov(data, tpi, 'time')
cov

cov['lags'] = lags
cov

fileout = os.path.join(os.getenv("SCRATCH"), 'covariance_OOPE_anomalies_tpi_sizes.nc')
fileout

delayed = cov.to_netcdf(fileout, compute=False)

with ProgressBar():
    delayed.compute()


