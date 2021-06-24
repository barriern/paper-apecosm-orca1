# # Step 2: Recovering detrended simulated SST anomalies
#
# ## Extracting the anomalies

import xarray as xr
import numpy as np
import apecosm.ts as ts
import scipy.signal as sig
import os.path

# **Note: here, the NEMO data has been extracted using `ncks -d olevel,0 ...`**

data = xr.open_mfdataset("data/nemo/*nc")
sst = data['thetao']

sst = np.squeeze(sst.to_masked_array())
print(sst.shape)
clim, sst = ts.get_monthly_clim(sst)  # time, lat, lon
del(clim) 
sst = sst.T  # lon, lat, time

# ## Detrending the time-series
#
# Here, SST time-series has been transposed in order to improve computation time.

nx, ny, nt = sst.shape

for i in range(nx):
    for j in range(ny):
        try:
            sst[i, j] = sig.detrend(sst[i, j])
        except:
            pass

# ## Writting the output dataset

dsout = xr.Dataset()
dsout['sst'] = (('x', 'y', 'time'), sst)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('modsst_anoms.nc')
