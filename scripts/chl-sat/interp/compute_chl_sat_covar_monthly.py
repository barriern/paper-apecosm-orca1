import xarray as xr
import scipy.signal as sig
import numpy as np
import os.path

import sys
sys.path.append("../../nino")
from extract_nino import read_index

daten, nino = read_index(filename='/home/barrier/Work/apecosm/ORCA1/DATA_APE_ORCA1/apecosm/index/oni.data', keepnan=False)

data = xr.open_mfdataset("data/anom*nc", combine='by_coords')
lat = data['lat'].values
lon = data['lon'].values
date = data['time.year'] * 100 + data['time.month']
date = date.values

dmax = np.min([date.max(), daten.max()])
dmin = np.max([date.min(), daten.min()])
dmax = 201812

inino = np.nonzero((daten <= dmax) & (daten >= dmin))
idata = np.nonzero((date <= dmax) & (date >= dmin))[0]

daten = daten[inino]
nino = nino[inino]
date = date[idata]
data = data.isel(time=idata)
data = data['chl_anom'].values

print(date)
print(daten)

data = data.T  # time, lat, lon -> lon, lat, time

nlon, nlat, ntime = data.shape

output = np.zeros((nlon, nlat))

for i in range(nlon):
    for j in range(nlat):
        temp = data[i, j]
        iok = np.nonzero(temp == temp)[0]
        if len(iok) > 0:
            output[i, j] = np.cov(temp[iok], nino[iok])[0, 1]

output = output.T  # lat, lon

dsout = xr.Dataset({'lon': lon, 'lat':lat})
dsout['cov'] = (['lat', 'lon'], output)
dsout.attrs['file'] = os.path.relpath(__file__)
dsout.to_netcdf("covariance_satellite_data.nc")





