# +
import xarray as xr
import numpy as np
import sys
import matplotlib.pyplot as plt
import os
sys.path.append(os.path.join('..', '..'))
from extract_nino import read_index

dirin = os.path.join('..', '..', 'data')
# -

dnino, nino = read_index(os.path.join(dirin, 'external', 'oni.data'))

iok = np.nonzero((dnino >= 195801) & (dnino <= 201812))
dnino = dnino[iok]
nino = nino[iok]

data = xr.open_dataset('hadsst_anoms.nc')
lon = data['lon']
lat = data['lat']
sst = data['sst'].to_masked_array().T
data

nx, ny, nt = sst.shape
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
dsout.to_netcdf('cov_hadsst_oni_tpi.nc')
