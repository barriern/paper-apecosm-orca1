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

data = xr.open_dataset('data/satellite_anoms_%s.nc' %varname)
data = data[varname]
data
# -

dates = data['time.year'].values * 100 + data['time.month'].values

ndate, nino = read_index('../../../data/index/oni.data')

datemin = np.max([ndate.min(), dates.min()])
datemin

datemax = np.min([ndate.max(), dates.max()])
datemax

inino = np.nonzero((ndate <= datemax) & (ndate >= datemin))[0]

isat = np.nonzero((dates <= datemax) & (dates >= datemin))[0]

data = data.isel(time=isat)
data

ninods = xr.Dataset()
ninods['time'] = data['time']
ninods['nino'] = (['time'], nino[inino])
ninods

cov = xr.cov(ninods['nino'], data, 'time')
cov.to_netcdf('data/sat_covariance_nino34_%s.nc' %varname)
