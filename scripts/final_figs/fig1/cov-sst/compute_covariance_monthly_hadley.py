import xarray as xr
import numpy as np
import sys
sys.path.append('../nino')
from extract_nino import read_index
import matplotlib.pyplot as plt
import os

dnino, nino = read_index('../data/index/oni.data')

iok = np.nonzero((dnino >= 195801) & (dnino <= 201812))
dnino = dnino[iok]
nino = nino[iok]

tpi = xr.open_dataset('../data/filt_tpi.nc')
dtpi = tpi['time'].values
tpi = tpi['tpi_filt'].values

itpi = np.nonzero(np.isnan(tpi) == False)
tpi = tpi[itpi]
#tpi = (tpi - tpi[:].mean()) / np.std(tpi[:])

data = xr.open_dataset('hadsst_anoms.nc')
lon = data['lon']
lat = data['lat']
sst = data['sst'].to_masked_array()

nx, ny, nt = sst.shape
covtpi = np.zeros((nx, ny))
covoni = np.zeros((nx, ny))

for i in range(nx):
    for j in range(ny):
        temp = sst[i, j]
        covtpi[i, j] = np.cov(temp[itpi], tpi)[0, 1]
        covoni[i, j] = np.cov(temp, nino)[0, 1]

dsout = xr.Dataset()
dsout['lon'] = lon
dsout['lat'] = lat
dsout['covtpi'] = (['lon', 'lat'], covtpi)
dsout['covoni'] = (['lon', 'lat'], covoni)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('cov_hadsst_oni_tpi.nc')

