import xarray as xr
import numpy as np
import apecosm.ts as ts
import scipy.signal as sig
import os.path

data = xr.open_mfdataset("data/CHL*nc")
chl = data['NCHL'] + data['DCHL']

chl = np.squeeze(chl.to_masked_array())
print(chl.shape)
clim, chl = ts.get_monthly_clim(chl)  # time, lat, lon
del(clim) 
chl = chl.T  # x,.

nx, ny, nt = chl.shape

for i in range(nx):
    for j in range(ny):
        try:
            chl[i, j] = sig.detrend(chl[i, j])
        except:
            pass

dsout = xr.Dataset()
dsout['chl'] = (('x', 'y', 'time'), chl)
#dsout['time'] = time
#dsout['lat'] = (['lat'], lat)
#dsout['lon'] = (['lon'], lon)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('modchl_anoms.nc')
