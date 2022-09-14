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
import numpy as np

varname = 'growthTrend'
# -

filelist = []
for y in [1982, 1997, 2015]:
    pattern = '/home1/scratch/nbarrier/*%s_Y%.4d*nc' %(varname, y)
    print(pattern)
    temp = glob(pattern)[0]
    filelist.append(temp)
filelist

data = xr.open_mfdataset(filelist)[varname]
data

clim = xr.open_dataset('data/pacific_clim_%s.nc' %varname)[varname]
clim

anom = data.groupby('time.month') - clim
anom

months = anom['time.month'].values
iok = np.nonzero(months >=10)[0]
iok

anom = anom.isel(time=iok)
anom['time']

compo = anom.mean(dim='time')
compo

fileout = 'data/composite_%s_map.nc' %varname
fileout


compo.to_netcdf(fileout)
