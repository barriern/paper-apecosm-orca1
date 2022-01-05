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
from glob import glob
import xarray as xr
import os
from dask.diagnostics import ProgressBar

grid = 'v_active'
datestart = '1971-01-01'
dateend = '2000-12-31'
# -

dirin = os.getenv('SCRATCH')
filelist = glob('%s/*JRA*%s*' %(dirin, grid))
filelist.sort()
filelist

data = xr.open_mfdataset(filelist)
data

data = data.sel(time=slice(datestart, dateend))
data

clim = data.groupby('time.month').mean(dim='time')
clim

clim.attrs['description'] = 'Monthly climatology'
clim.attrs['datestart'] = datestart
clim.attrs['dateend'] = dateend
clim

outfile = '%s/pacific_clim_%s.nc' %(dirin, grid)
outfile

delayed = clim.to_netcdf(outfile, compute=False)
with ProgressBar():
    delayed.compute()


