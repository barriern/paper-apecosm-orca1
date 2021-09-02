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

grid = 'mdiff_trend'
datestart = '1996-01-01'
dateend = '2001-12-31'
# -

dirin = os.getenv('SCRATCH')
filelist = glob('%s/*JRA*%s*' %(dirin, grid))
filelist.sort()
filelist

data = xr.open_mfdataset(filelist)
data

data = data.sel(time=slice(datestart, dateend))
data

data.attrs['description'] = 'Monthly dataatology'
data.attrs['datestart'] = datestart
data.attrs['dateend'] = dateend
data

outfile = '%s/pacific_nino97_%s.nc' %(dirin, grid)
outfile

delayed = data.to_netcdf(outfile, compute=False)
with ProgressBar():
    delayed.compute()


