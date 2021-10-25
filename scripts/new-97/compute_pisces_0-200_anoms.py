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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar

plt.rcParams['text.usetex'] = False

zmax = 50
varname = 'ZOO2'

datestart = '1997-01'
dateend = '1999-12'
# -

# ## Extracting the climatology
#

clim = xr.open_dataset('data/pacific_0-%d_%s_clim.nc' %(zmax, varname))
clim

varclim = clim[varname]
varclim

# ## Extracting the data

data = xr.open_dataset('data/pacific_0-%d_%s.nc' %(zmax, varname)).sel(time_counter=slice(datestart, dateend))
data

var = data[varname]
var

# ## Computing the anomalies

anom = var.groupby('time_counter.month') - varclim
anom

anom.to_netcdf('data/pacific_0-%d_%s_anom.nc' %(zmax, varname))
