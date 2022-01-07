# # Step 1/4: Processing Hadley SST
#
# **Objectives: extract the detrended SST anomalies from the Hadley dataset over the 1958-2018 period**

import xarray as xr
import numpy as np
import apecosm.ts as ts
import scipy.signal as sig
import os.path

data = xr.open_mfdataset("/home1/scratch/nbarrier/HadISST_sst.nc")
data

# ## Recovering the anomalies

# Recovering the date of the Hadley dataset and extracts the overlapping period (1958-2018)

year = data['time.year'].values
month = data['time.month'].values
date = year * 100 + month 
print(date)

iok = np.nonzero((date >= 195801) & (date <= 201812))[0]

lon = data['longitude'].values
lat = data['latitude'].values

data = data.isel(time=iok)
sst = data['sst']
dims = sst.dims  # time, lat, lon
print(sst.shape)
time = data['time']
sst

clim = data.sel(time=slice('1971-01-01','2000-12-31'))
clim = clim['sst'].groupby('time.month').mean(dim='time')
clim

anoms = data['sst'].groupby('time.month') - clim
anoms

# ## Detrending the anomalies

sst = anoms.to_masked_array().T
sst.shape

nx, ny, nt = sst.shape

for i in range(nx):
    for j in range(ny):
        try:
            sst[i, j] = sig.detrend(sst[i, j])
        except:
            pass

# ## Writting the output

dsout = xr.Dataset()
dsout['sst'] = (('time', 'lat', 'lon'), sst.T)
dsout['time'] = time
dsout['lat'] = (['lat'], lat)
dsout['lon'] = (['lon'], lon)
dsout.to_netcdf('../data/hadsst_anoms.nc')
