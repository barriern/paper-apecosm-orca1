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

varname = 'uo'
nino_years = [1982, 1997, 2015]
nino_dates = []
for y in nino_years:
    for m in range(10, 13):
        nino_dates.append(y * 100 + m)
nino_dates = np.array(nino_dates)
# -

data = xr.open_dataset('data/equatorial_full_%s.nc' %varname)[varname]
data

clim = data.sel(time_counter=slice('1971-01-01', '2000-12-31')).groupby('time_counter.month').mean()
clim

anom = data.groupby('time_counter.month') - clim
anom

dates = data['time_counter.year'].values * 100 + data['time_counter.month'].values 
test = np.array([d in nino_dates for d in dates])
iok = np.nonzero(test)[0]
dates[iok]

anom = anom.isel(time_counter=iok)
anom['time_counter']

compo = anom.mean(dim='time_counter')
compo

compo.to_netcdf('data/composite_%s_profile.nc' %varname)
compo.plot()

mean = clim.mean(dim='month')
mean.to_netcdf('data/mean_%s_profile.nc' %varname)


