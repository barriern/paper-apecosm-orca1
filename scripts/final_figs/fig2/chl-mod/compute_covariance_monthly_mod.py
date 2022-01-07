import xarray as xr
import numpy as np
import sys
sys.path.append('../../../nino')
from extract_nino import read_index
import matplotlib.pyplot as plt
import os
import scipy.signal as sig

dnino, nino = read_index('../../../data/index/oni.data')

iok = np.nonzero((dnino >= 195801) & (dnino <= 201812))
dnino = dnino[iok]
nino = nino[iok]

data = xr.open_dataset('data/anom_chl_monthly_model.nc')
#don = data['lon']
#lat = data['lat']
chl = data['anom_chl'].to_masked_array()
chl[chl.mask] = 0
chl = sig.detrend(chl.T)
chl.shape

nx, ny, nt = chl.shape
covtpi = np.zeros((nx, ny))
covoni = np.zeros((nx, ny))

for i in range(nx):
    for j in range(ny):
        temp = chl[i, j]
        covoni[i, j] = np.cov(temp, nino)[0, 1]

dsout = xr.Dataset()
#dsout['lon'] = lon
#dsout['lat'] = lat
dsout['covoni'] = (['y', 'x'], covoni.T)
dsout.to_netcdf('data/cov_modchl_oni_tpi.nc')



