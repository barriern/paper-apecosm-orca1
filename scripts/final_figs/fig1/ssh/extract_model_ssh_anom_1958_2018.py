# # Step 2: Recovering detrended simulated SST anomalies
#
# ## Extracting the anomalies

import xarray as xr
import numpy as np
import apecosm.ts as ts
import scipy.signal as sig
import os.path

# **Note: here, the NEMO data has been extracted using `ncks -d olevel,0 ...`**

data = xr.open_mfdataset("/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/JRA_CO2/*ssh*nc")
sst = data['zos']
sst

clim = sst.sel(time_counter=slice('1971-01-01', '2000-12-31'))
clim = clim.groupby('time_counter.month').mean(dim='time_counter')
clim

anom = sst.groupby('time_counter.month') - clim
anom

# ## Detrending the time-series
#
# Here, SST time-series has been transposed in order to improve computation time.

sst = anom.to_masked_array().T

nx, ny, nt = sst.shape
sst.shape

for i in range(nx):
    for j in range(ny):
        try:
            sst[i, j] = sig.detrend(sst[i, j])
        except:
            pass

# ## Writting the output dataset

dsout = xr.Dataset()
dsout['ssh'] = (('x', 'y', 'time'), sst)
dsout.to_netcdf('../data/modssh_anoms.nc')
