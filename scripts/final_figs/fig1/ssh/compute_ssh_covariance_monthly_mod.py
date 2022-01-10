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

data = xr.open_dataset('../data/modssh_anoms.nc')
#don = data['lon']
#lat = data['lat']
sst = data['zos']
sst = sst.sel(time_counter=slice('1993-01-01', '2018-12-31'))

ny, nx, nt = sst.shape
covoni = np.zeros((ny, nx))
sst.shape

for j in range(ny):
    for i in range(nx):
        temp = sst.isel(y=j, x=i).values
        covoni[j, i] = np.cov(temp, nino)[0, 1]

dsout = xr.Dataset()
dsout['covoni'] = (['y', 'x'], covoni)
dsout.to_netcdf('../data/cov_modssh_oni_tpi.nc')

