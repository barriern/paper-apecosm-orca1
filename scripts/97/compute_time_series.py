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

# # Time-serie computation
#
# Time-series are computed for the central Pacific (entre 200 et 250) and the western Pacific (150 et 180), and between 5N and 5S
#
# ## Extraction of the domain mask

# +
import xarray as xr

latmax = 5
varname = 'madv_trend'

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/processed_pacific'

mesh = xr.open_dataset('%s/pacific_mesh_mask.nc' %dirin).isel(z=0)
mesh
# -

tmask = mesh['tmaskutil']
tmask

surf = mesh['e1t'] * mesh['e2t'] * tmask
surf

lat = mesh['gphit']
lon = mesh['glamt']
lon

lon = (lon + 360) % 360

mask_central = ((abs(lat) <= latmax) & (lon >= 200) & (lon <= 250)).astype(int)
(mask_central + tmask).plot()

mask_west = ((abs(lat) <= latmax) & (lon >= 150) & (lon <= 180)).astype(int)
(mask_west + tmask).plot()

# ## Loading the data and computing the mean anomalies

dataclim = xr.open_dataset('%s/pacific_clim_%s.nc' %(dirin, varname))
dataclim

varclim = dataclim[varname]
varclim

data = xr.open_dataset('%s/pacific_nino97_%s.nc' %(dirin, varname))
data

var = data[varname]
var

# Now, we compute the anomalies:

var = var.groupby('time.month') - varclim

ts_central = (var * surf * mask_central).sum(dim=['x', 'y']) / (surf * mask_central).sum(dim=['y', 'x'])
ts_central.name = '%s_central' %varname
ts_central

ts_west = (var * surf * mask_west).sum(dim=['x', 'y']) / (surf * mask_west).sum(dim=['y', 'x'])
ts_west.name = '%s_west' %varname
ts_west

dataout = xr.Dataset()
dataout['%s_west' %varname] = ts_west
dataout['%s_central' %varname] = ts_central
dataout

fileout = '%s/%s_timeseries.nc' %(dirin, varname)
fileout

dataout.to_netcdf(fileout)


