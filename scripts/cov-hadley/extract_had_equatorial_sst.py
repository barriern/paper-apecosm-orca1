import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import os

data = xr.open_dataset('data/HadISST_sst.nc')
year = data['time.year'].values
itime = np.nonzero((year >= 1958) & (year <= 2018))[0]
data = data.isel(time=itime)


lon = data['longitude'].values
lat = data['latitude'].values
nlat = len(lat)
nlon = len(lon)


lon2d, lat2d = np.meshgrid(lon, lat)
weights = np.cos(np.deg2rad(lat2d))
weights = weights[np.newaxis, :, :]

lonmin = 150
lonmax = -80
testlat = (np.abs(lat2d) <= 2)
testlon = (lon2d <= lonmax) | (lon2d >= lonmin)
test = testlat & testlon

ilat, ilon = np.nonzero(test)

output = np.zeros((nlat, nlon))
output[ilat, ilon] = 1

sst = data['sst'].to_masked_array()
mask = 1 - np.ma.getmaskarray(sst)  # mask = 1 si masque, ici mask = 1 si non masque

den = (weights * mask)[:, ilat, ilon].sum(axis=-1)
num = (weights * mask * sst)[:, ilat, ilon].sum(axis=-1)

output = num / den
print(output.shape)

dsout = xr.Dataset()
dsout['time'] = data['time']
dsout['sst'] = (['time'], output)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('obs_equatorial_mean.nc')
