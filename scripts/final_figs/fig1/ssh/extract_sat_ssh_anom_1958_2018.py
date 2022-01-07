# # Step 2: Recovering detrended simulated SST anomalies
#
# ## Extracting the anomalies

import xarray as xr
import numpy as np
import scipy.signal as sig
import os.path

# **Note: here, the NEMO data has been extracted using `ncks -d olevel,0 ...`**

data = xr.open_mfdataset("../data/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M-m_1641569308380.nc")
sst = data['sla']

clim = sst.sel(time=slice('1993-01-01', '2018-12-31'))
clim = clim.groupby('time.month').mean(dim='time')
clim

anom = sst.groupby('time.month') - clim
anom

# ## Detrending the time-series
#
# Here, SST time-series has been transposed in order to improve computation time.

sst = anom.to_masked_array().T
sst.shape

nx, ny, nt = sst.shape

for i in range(nx):
    for j in range(ny):
        try:
            sst[i, j] = sig.detrend(sst[i, j])
        except:
            pass

# ## Writting the output dataset

dsout = xr.Dataset()
dsout['ssh'] = (('x', 'y', 'time'), sst)
dsout['time'] = anom['time']
dsout.to_netcdf('../data/obsssh_anoms.nc')
