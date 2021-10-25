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

# # Computation of full time-series anomalies
#
# This will be used for computation of covariances and EOF analysis

# +
import xarray as xr
from glob import glob
import os

varname = 'OOPE'

dirin = os.getenv('SCRATCH')
dirin
# -

# ## Reading the climatology

clim = xr.open_dataset('data/pacific_clim_%s.nc' %varname)
clim

varclim = clim[varname]
varclim

# ## Reading the entire time-series to extract anomalies

filepattern = os.path.join(dirin, 'pacific*%s*nc' %varname)
filepattern

filelist = glob(filepattern)
filelist.sort()
filelist

for f in filelist:
    
    print('Processing file ', f)
    data = xr.open_dataset(f)
    data

    var = data[varname]
    var

    anom = var.groupby('time.month') - varclim
    anom.name = varname
    anom

    fileout = f.replace('pacific', 'anomalies_pacific')
    fileout

    anom.to_netcdf(fileout, unlimited_dims='time')


