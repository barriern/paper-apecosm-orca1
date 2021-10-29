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

dirin = os.getenv('DATAWORK')

latmax = 30
lonmin = 150
lonmax = -120

dirin = os.path.join(dirin, 'apecosm/apecosm_orca1/processed_pacific/')
# -

# ## Loading EOF Pcs

fileout = '%s/full_eof_pacific_OOPE_latmax_%d_lonmin_%d_lonmax_%d.nc' %(dirin, latmax, lonmin, lonmax)
fileout

eof = xr.open_dataset(fileout)
eof

pcs = eof['eofpcs']
pcs.load()
pcs

# ## Loading the datasets

dirin = os.getenv('SCRATCH')
data = xr.open_dataset(os.path.join(dirin, 'detrended_pacific_OOPE_anomalies.nc'))
data = data['OOPE']
data  # y, x, w, time

ny, nx, nw, ntime = data.shape


# ## Preparing stuff for computation of lead-lag covariances

def gufunc_cov(x, y):
    print('X.shape ', x.shape)  # (6, 40, 100, 1, 732)  => (x, y, w, 1, ntime)
    print('Y.shape ', y.shape)  #(100, 2, 732) => (weight, eof, time)
    output = (x - np.mean(x, axis=-1, keepdims=True)) * (y - np.mean(y, axis=-1, keepdims=True)).mean(dim=-1)
    output /= (x.shape[-1] - 1)

    return output


def xarray_cov(x, y, dim):
    return xr.apply_ufunc(
        gufunc_cov,
        x,
        y,
        input_core_dims=[[dim], [dim]],
        #output_core_dims=[['lags']],
        dask="parallelized",
        output_dtypes=[np.float32],
        #dask_gufunc_kwargs = {'output_sizes' : {'lags': nlags}}
    )


data = data.chunk({'w' : 20, 'x' : 30, 'y' : 40})
data

cov = xr.cov(data, pcs, 'time')
cov

fileout = os.path.join(os.getenv("SCRATCH"), 'covariance_OOPE_anomalies_pcs_latmax_%d_lonmin_%d_lonmax_%d_all_sizes.nc' %(latmax, lonmin, lonmax))
fileout

delayed = cov.to_netcdf(fileout, compute=False)

with ProgressBar():
    delayed.compute()
