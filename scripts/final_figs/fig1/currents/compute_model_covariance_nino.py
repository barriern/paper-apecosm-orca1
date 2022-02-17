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
import sys
sys.path.append('../../../nino/')
from extract_nino import read_index
import numpy as np

varname = 'uo'

data = xr.open_dataset('data/model_anoms_%s.nc' %varname)
data = data[varname]
data
# -

dates = data['time_counter.year'].values * 100 + data['time_counter.month'].values

ndate, nino = read_index('../../../data/index/oni.data')

datemin = np.max([ndate.min(), dates.min()])
datemin

datemax = np.min([ndate.max(), dates.max()])
datemax

inino = np.nonzero((ndate <= datemax) & (ndate >= datemin))[0]

isat = np.nonzero((dates <= datemax) & (dates >= datemin))[0]

data = data.isel(time_counter=isat)
data

ninods = xr.Dataset()
ninods['time_counter'] = data['time_counter']
ninods['nino'] = (['time_counter'], nino[inino])
ninods

cov = xr.cov(ninods['nino'], data, 'time_counter')
cov.to_netcdf('data/model_covariance_nino34_%s.nc' %varname)


