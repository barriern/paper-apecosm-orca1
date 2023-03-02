# # Step 2: Recovering detrended simulated SST anomalies
#
# ## Extracting the anomalies

# +
import xarray as xr
import numpy as np
import scipy.signal as sig
import os.path
from dask.diagnostics import ProgressBar

dirin = os.path.join('..', '..', 'data', 'external')
# -

# **Note: here, the NEMO data has been extracted using `ncks -d olevel,0 ...`**

data = xr.open_dataset(os.path.join(dirin, 'cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M-m_1641569308380.nc'))
sst = data['sla']
sst = sst.chunk({'latitude': 100, 'longitude': 100})
sst

clim = sst.sel(time=slice('1993-01-01', '2018-12-31'))
clim = clim.groupby('time.month').mean(dim='time')
clim

anom = sst.groupby('time.month') - clim
anom = anom.chunk({'time' : anom.shape[0]})
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

detrended = xarray_detrend(anom, dim='time')
detrended = detrended.compute()

# ## Writting the output dataset

detrended.to_netcdf('obsssh_anoms.nc')
