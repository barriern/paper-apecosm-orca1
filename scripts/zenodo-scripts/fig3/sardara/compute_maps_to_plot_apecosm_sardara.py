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

import xarray as xr
import numpy as np

# +
dates_nino = []
for y in [2009, 2014, 2015]:
    for m in [10, 11, 12]:
        dates_nino.append(y * 100 + m)
    for m in [1, 2, 3]:
        dates_nino.append((y + 1) * 100 + m)
    
dates_nino
# -

dates_nina = []
for y in [2007, 2008, 2010, 2011]:
    for m in [10, 11, 12]:
        dates_nina.append(y * 100 + m)
    for m in [1, 2, 3]:
        dates_nina.append((y + 1) * 100 + m)
dates_nina

data1 = xr.open_mfdataset('sardara/regridded_catch_gear_PS_species_SKJ_1x1.nc')
data2 = xr.open_mfdataset('sardara/regridded_catch_gear_PS_species_YFT_1x1.nc')
data = data1['catch'] + data2['catch']
data

dates_ape = data['time'].values
dates_ape[:5]

inino = np.nonzero(np.array([d in dates_nino for d in dates_ape]))[0]
dates_ape[inino].shape

inina = np.nonzero(np.array([d in dates_nina for d in dates_ape]))[0]
dates_ape[inina]

compo_nino = data.isel(time=inino).mean(dim='time')
compo_nino

compo_nina = data.isel(time=inina).mean(dim='time')
compo_nina

compo_nina.to_netcdf('compo_sardara_nina.nc')

compo_nino.to_netcdf('compo_sardara_nino.nc')
