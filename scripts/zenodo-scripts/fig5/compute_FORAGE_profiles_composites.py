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

varname = 'FORAGE'
nino_years = [1982, 1997, 2015]
nino_dates = []
for y in nino_years:
    for m in range(10, 13):
        nino_dates.append(y * 100 + m)
nino_dates = np.array(nino_dates)
# -

data = xr.open_dataset('../fig4/equatorial_full_%s.nc' %varname)[varname]
data

clim = data.sel(time=slice('1971-01-01', '2000-12-31')).groupby('time.month').mean()
clim

anom = data.groupby('time.month') - clim
anom

dates = data['time.year'].values * 100 + data['time.month'].values 
test = np.array([d in nino_dates for d in dates])
iok = np.nonzero(test)[0]
dates[iok]

anom = anom.isel(time=iok)
anom['time']

compo = anom.mean(dim='time')
compo

compo.to_netcdf('composite_%s_profile.nc' %varname)

mean = clim.mean(dim='month')
mean.to_netcdf('mean_%s_profile.nc' %varname)


