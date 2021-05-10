import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import apecosm.ts as ts
import os.path

mesh = xr.open_dataset('../../data/mesh_mask_eORCA1_v2.2.nc')
tmask = mesh['tmask'].values[0][0]
lon = mesh['glamt'].values[0]
lat = mesh['gphit'].values[0]
e1t = mesh['e1t'].values
e2t = mesh['e2t'].values

output = tmask.copy()

lonmin = 150
lonmax = -80

test = (np.abs(lat) <= 2)
test2 = (lon <= lonmax) | (lon >= lonmin)
test = (test & test2 & tmask)

ilat, ilon = np.nonzero(test == True)
output[ilat, ilon] = 2

surf = e1t * e2t

data = xr.open_mfdataset('data/CHL*nc')
time = data['time_counter']
data = data['NCHL'] + data['DCHL']
data = np.squeeze(data.to_masked_array())

num = (surf * data)[:, ilat, ilon]
den = (surf)[:, ilat, ilon]

output = np.sum(num, axis=-1) / np.sum(den, axis=-1)

dsout = xr.Dataset()
dsout['time'] = time
dsout['chl'] = (['time_counter'], output)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('simulated_equatorial_mean.nc')
