import xarray as xr
import numpy as np
import sys
sys.path.append('../../../nino')
from extract_nino import read_index
import matplotlib.pyplot as plt
import os

dnino, nino = read_index('../../../data/index/oni.data')

iok = np.nonzero((dnino >= 195801) & (dnino <= 201812))
dnino = dnino[iok]
nino = nino[iok]

data = xr.open_dataset('../data/hadsst_anoms.nc')
lon = data['lon']
lat = data['lat']
sst = data['sst'].to_masked_array().T
data

nx, ny, nt = sst.shape
covtpi = np.zeros((nx, ny))
covoni = np.zeros((nx, ny))
sst.shape

for i in range(nx):
    for j in range(ny):
        temp = sst[i, j]
        covoni[i, j] = np.cov(temp, nino)[0, 1]

dsout = xr.Dataset()
dsout['lon'] = lon
dsout['lat'] = lat
dsout['covoni'] = (['lon', 'lat'], covoni)
dsout.to_netcdf('../data/cov_hadsst_oni_tpi.nc')



