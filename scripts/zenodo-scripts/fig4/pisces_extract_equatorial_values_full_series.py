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
from glob import glob
from dask.diagnostics import ProgressBar

dirin = '../data/nemo-pisces'
grid = 'speed_U'
varname = 'uo'
ilat = 61
# -

mesh = xr.open_dataset('../data/static/pacific_mesh_mask.nc').isel(y=ilat, z=0)
lon0 = mesh['glamt']
lon0

filelist = glob('%s/pacific_nico*%s*nc' %(dirin, grid))
filelist.sort()
filelist

data = xr.open_mfdataset(filelist, decode_times=False)[varname].isel(y=ilat)
data

data['x'] = (lon0 + 360) % 360
data = data.sortby(data['x'])
data

delayed = data.to_netcdf('equatorial_full_%s.nc' %varname, compute=False)

with ProgressBar():
    delayed.compute()


