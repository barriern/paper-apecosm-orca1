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

varname = 'uo'
grid = 'speed_U'
write_meshout = True

import xarray as xr
import os
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar, visualize, Profiler, ResourceProfiler, CacheProfiler
plt.rcParams['text.usetex'] = False

# ## Extracting the domain

mesh = xr.open_dataset('../data/mesh_mask_eORCA1_v2.2.nc').isel(t=0)
mesh

lon = mesh['glamt']
lon

lat = mesh['gphit']
lat

ilat = (abs(lat) <= 40)
ilat

ilon = (lon >= 120) | (lon <= -60)
ilon

domain = ilon & ilat
domain

domain.plot()

meshout = mesh.where(domain == True, drop=True)
meshout

if(write_meshout):
    print('write mesh out')
    meshout.to_netcdf('extracted_mesh_file.nc')

# ## Extracting the data

data = xr.open_mfdataset('data/surface*%s*nc' %grid).isel(olevel=0)
data

var = data[varname]
var

var = var.where(domain==True, drop=True)
var

# ## Extracting the climatology

clim = xr.open_dataset('clim_%s.nc' %varname)
clim

varclim = clim[varname]
varclim

varclim = varclim.where(domain==True, drop=True)
varclim

anom = var - varclim
anom

delayed = anom.to_netcdf('anom_%s.nc' %varname, compute=False)
delayed

with ProgressBar():
    delayed.compute()


