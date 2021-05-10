import xarray as xr
import numpy as np
import apecosm.ts as ts
import scipy.signal as sig
import os.path

data = xr.open_mfdataset("data/nemo/*nc")
sst = data['thetao']

sst = np.squeeze(sst.to_masked_array())
print(sst.shape)
clim, sst = ts.get_monthly_clim(sst)  # time, lat, lon
del(clim) 
sst = sst.T  # x,.

nx, ny, nt = sst.shape

for i in range(nx):
    for j in range(ny):
        try:
            sst[i, j] = sig.detrend(sst[i, j])
        except:
            pass

dsout = xr.Dataset()
dsout['sst'] = (('x', 'y', 'time'), sst)
#dsout['time'] = time
#dsout['lat'] = (['lat'], lat)
#dsout['lon'] = (['lon'], lon)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('modsst_anoms.nc')
