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
lon = data['longitude']
lat = data['latitude']
sst = data['sla']
sst = sst.sel(time=slice('1993-01-01', '2018-12-31'))
sst

ny, nx, nt = sst.shape
covoni = np.zeros((ny, nx))

for j in range(ny):
    for i in range(nx):
        temp = sst.isel(longitude=i, latitude=j).values
        covoni[j, i] = np.cov(temp, nino)[0, 1]

dsout = xr.Dataset()
dsout['lon'] = lon
dsout['lat'] = lat
dsout['covoni'] = (['y', 'x'], covoni)
dsout.to_netcdf('../data/cov_obsssh_oni_tpi.nc')

