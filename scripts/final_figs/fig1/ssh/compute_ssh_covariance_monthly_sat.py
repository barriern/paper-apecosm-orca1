import xarray as xr
import numpy as np
import sys
sys.path.append('../../../nino')
from extract_nino import read_index
import matplotlib.pyplot as plt
import os

dnino, nino = read_index('../../../data/index/oni.data')

iok = np.nonzero((dnino >= 199301) & (dnino <= 201812))
dnino = dnino[iok]
nino = nino[iok]
nino.shape

data = xr.open_dataset('../data/obsssh_anoms.nc')
#don = data['lon']
#lat = data['lat']
sst = data['ssh']
sst
sst = sst.sel(time=slice('1993-01-01', '2018-12-31'))
sst

nx, ny, nt = sst.shape
covoni = np.zeros((nx, ny))
sst.shape

for i in range(nx):
    for j in range(ny):
        temp = sst[i, j]
        covoni[i, j] = np.cov(temp, nino)[0, 1]

dsout = xr.Dataset()
#dsout['lon'] = lon
#dsout['lat'] = lat
dsout['covoni'] = (['x', 'y'], covoni)
dsout.to_netcdf('../data/cov_obsssh_oni_tpi.nc')

