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

varname = 'OOPE'
# -

const = xr.open_dataset('../data/static/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100
const

filelist = []
for y in [1982, 1997, 2015]:
    pattern = '../data/apecosm/*%s_Y*%d*' %(varname, y)
    print(pattern)
    temp = glob(pattern)[0]
    filelist.append(temp)
filelist

data = xr.open_mfdataset(filelist)[varname]
data = data.rename({'w': 'l'})
data['l'] = const['l']
data = data.sel(l=[3, 20, 90], method='nearest')
data

clim = xr.open_dataset('clim_epi_%s.nc' %varname)[varname]
clim

anom = data.groupby('time.month') - clim
anom

# +
months = anom['time.month'].values
iok = np.nonzero(months >=10)[0]
print(iok)

anom = anom.isel(time=iok)
print(anom['time'])

compo = anom.mean(dim='time')
print(compo)

fileout = 'composite_%s_map.nc' %varname
print(fileout)

compo.to_netcdf(fileout)
# -


