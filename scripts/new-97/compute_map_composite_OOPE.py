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
import xarray as xr
from glob import glob
import numpy as np

filelist = []
for y in [1982, 1997, 2015]:
    temp = glob('/home1/scratch/nbarrier/*OOPE_Y%.4d*nc' %y)[0]
    filelist.append(temp)
# -

data = xr.open_mfdataset(filelist)['OOPE']
data

clim = xr.open_dataset('data/pacific_clim_OOPE.nc')['OOPE']
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

compo.to_netcdf('data/composite_OOPE_map.nc')
