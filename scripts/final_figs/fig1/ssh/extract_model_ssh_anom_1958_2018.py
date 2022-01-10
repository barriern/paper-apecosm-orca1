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
sst = sst.chunk({'x': 100, 'y': 100})

clim = sst.sel(time_counter=slice('1993-01-01', '2018-12-31'))
clim = clim.groupby('time_counter.month').mean(dim='time_counter')
clim

anom = sst.groupby('time_counter.month') - clim
anom = anom.chunk({'time_counter': anom.shape[0]})
anom


# ## Detrending the time-series
#
# Here, SST time-series has been transposed in order to improve computation time.

# +
def gufunc_detrend(x):
    x[np.isnan(x)] = 0
    return sig.detrend(x)

def xarray_detrend(x, dim):
    return xr.apply_ufunc(
        gufunc_detrend,
        x,
        input_core_dims=[[dim]],
        output_core_dims=[[dim]],
        dask="parallelized",
        output_dtypes=[np.float32],
    )


# -

detrended = xarray_detrend(anom, dim='time_counter')
detrended

# ## Writting the output dataset

detrended.to_netcdf('../data/modssh_anoms.nc')
