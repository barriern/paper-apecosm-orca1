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
import scipy.signal as sig 
import numpy as np
from dask.diagnostics import ProgressBar
from glob import glob
import sys
import pandas as pd
import matplotlib.pyplot as plt

dirin = '../data/apecosm/'

latmax = 10
lonmin = 150
lonmax = -100
# -

# ## Loading ONI index

sys.path.append('../')
from extract_nino import read_index
dateoni, oni = read_index('../data/external/oni.ascii.txt')
len(dateoni)
plt.plot(oni)
oni = xr.DataArray(data=oni, dims='time')
oni

# ## Loading the datasets

dirin = os.getenv('SCRATCH')
data = xr.open_dataset(os.path.join(dirin, 'detrended_pacific_OOPE_anomalies.nc'))
data = data['OOPE']
data  # y, x, w, time

ny, nx, nw, ntime = data.shape

# ## Preparing stuff for computation of lead-lag covariances

data = data.chunk({'x' : 30, 'y' : 40})
data

data.shape, oni.shape

cov = xr.cov(data, oni, 'time')
cov.name = 'OOPE'

# +
# cov['lags'] = lags
# cov
# -

fileout = os.path.join(os.getenv("SCRATCH"), 'covariance_OOPE_anomalies_oni_sizes.nc')
fileout

delayed = cov.to_netcdf(fileout, compute=False)

with ProgressBar():
    delayed.compute()


