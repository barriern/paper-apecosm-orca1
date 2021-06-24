# # Step 1/4: Processing Hadley SST
#
# **Objectives: extract the detrended SST anomalies from the Hadley dataset over the 1958-2018 period**

import xarray as xr
import numpy as np
import apecosm.ts as ts
import scipy.signal as sig
import os.path

data = xr.open_mfdataset("data/HadISST_sst.nc")

# ## Recovering the anomalies

# Recovering the date of the Hadley dataset and extracts the overlapping period (1958-2018)

year = data['time.year'].values
month = data['time.month'].values
date = year * 100 + month 
print(date)

iok = np.nonzero((date >= 195801) & (date <= 201812))[0]
print(date[iok])
print(iok)

lon = data['longitude'].values
lat = data['latitude'].values

data = data.isel(time=iok)
sst = data['sst']
dims = sst.dims  # time, lat, lon
print(sst.shape)
time = data['time']

sst = sst.to_masked_array()
clim, sst = ts.get_monthly_clim(sst)  # time, lat, lon
del(clim) 
sst = sst.T  # lon, lat, time

# ## Detrending the anomalies

nx, ny, nt = sst.shape

for i in range(nx):
    for j in range(ny):
        try:
            sst[i, j] = sig.detrend(sst[i, j])
        except:
            pass

# ## Writting the output

dsout = xr.Dataset()
dsout['sst'] = (('lon', 'lat', 'time'), sst)
dsout['time'] = time
dsout['lat'] = (['lat'], lat)
dsout['lon'] = (['lon'], lon)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('hadsst_anoms.nc')
