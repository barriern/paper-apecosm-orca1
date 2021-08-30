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

# +
import xarray as xr
import matplotlib.pyplot as plt

grid = varname = 'OOPE'
latmax = 1
# -

# ## Reading mesh mask

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0)
mesh

lat = mesh['gphit']
lat

mesh = mesh.where(abs(lat) <= latmax)

lon = mesh['glamt'].mean(dim='y')
lon

lon = (lon + 360) % 360
lon

surf = mesh['e1t'] * mesh['e2t'] * mesh['tmask']
surf

# ## Loading the climatology

clim = xr.open_dataset('data/pacific_clim_%s.nc' %(grid))
clim

varclim = clim[varname]
varclim

# ## Loading the field

data = xr.open_dataset('data/pacific_nino97_%s.nc' %(grid))
data

var = data[varname]
var

# ## Computing the anomalies

anom = var.groupby('time.month') - varclim
anom

# ## Computing the meridional anomalies

anomout = (anom * surf).sum(dim='y') / surf.sum(dim='y')
anomout

anomout['x'] = lon

const  = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
lengths = const['length']
lengths

anomout['w'] = lengths.values * 100
anomout


