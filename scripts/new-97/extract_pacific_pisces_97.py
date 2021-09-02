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

grid = 'ptrc_T'
datestart = '1996-01-01'
dateend = '2001-12-31'
# -

dirin = os.getenv('SCRATCH')
filelist = glob('%s/*JRA*%s*' %(dirin, grid))
filelist.sort()
filelist

data = xr.open_mfdataset(filelist, decode_times=False)
data

if((grid == 'speed_U') | (grid == 'speed_V')):
    data = data.drop(['time_centered'])
else:       
    data = data.drop(['time_centered', 'time_counter_bounds'])
data

data = xr.decode_cf(data)
data['time_counter']

data = data.sel(time_counter=slice(datestart, dateend))
data['time_counter']

data.attrs['description'] = 'Extraction of specific dates'
data.attrs['datestart'] = datestart
data.attrs['dateend'] = dateend
data

outfile = '%s/pacific_nino97_%s.nc' %(dirin, grid)
outfile

delayed = data.to_netcdf(outfile, compute=False)
with ProgressBar():
    delayed.compute()


