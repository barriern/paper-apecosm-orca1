# # Step 2: Recovering detrended simulated SST anomalies
#
# ## Extracting the anomalies

# +
import xarray as xr
import numpy as np
import scipy.signal as sig
import os.path
from glob import glob

dirin = os.path.join('..', '..', 'data', 'nemo-pisces')
# -

# **Note: here, the NEMO data has been extracted using `ncks -d olevel,0 ...`**

filelist = glob(os.path.join(dirin, '*ssh*'))
filelist.sort()
filelist

data = xr.open_mfdataset(filelist)
sst = data['zos']

clim = sst.sel(time_counter=slice('1993-01-01', '2018-12-31'))
clim = clim.groupby('time_counter.month').mean(dim='time_counter')
clim

anom = sst.groupby('time_counter.month') - clim
anom = anom.chunk({'time_counter': anom.shape[0]})
anom


# ## Detrending the time-series

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

detrended.to_netcdf('modssh_anoms.nc')


