# +
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import os

dirin = ''
# -

data = xr.open_mfdataset('/home1/scratch/nbarrier/sat-chl/interpolated*nc')
year = data['time.year'].values

lon = data['lon'].values
lat = data['lat'].values
nlat = len(lat)
nlon = len(lon)

lon2d, lat2d = np.meshgrid(lon, lat)
weights = np.cos(np.deg2rad(lat2d))
weights = weights[np.newaxis, :, :]

lonmin = 150
lonmax = -80 + 360
testlat = (np.abs(lat2d) <= 2)
testlon = (lon2d <= lonmax) & (lon2d >= lonmin)
test = testlat & testlon


ilat, ilon = np.nonzero(test)

output = np.zeros((nlat, nlon))
output[ilat, ilon] = 1

chlor_a = data['chlor_a'].isel(time=0).values
count_a = data['count_a'].isel(time=0).values

chlor_a = np.ma.masked_where(count_a == 0, chlor_a)

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree()
ax = plt.gca(projection=proj)
cs = ax.pcolormesh(lon, lat, np.log10(chlor_a), transform=proj2)
ax.plot(lon2d[ilat, ilon], lat2d[ilat, ilon], marker='.', linestyle='none', transform=proj2)

chl = data['chlor_a'].to_masked_array()
count = data['count_a'].to_masked_array()

mask = (count >= 100 / 3.).astype(int)

den = (weights * mask)[:, ilat, ilon].sum(axis=-1)
num = (weights * mask * chl)[:, ilat, ilon].sum(axis=-1)

output = num / den
print(output.shape)

dsout = xr.Dataset()
dsout['time'] = data['time']
dsout['chl'] = (['time'], output)
dsout.to_netcdf('/home1/datawork/nbarrier/apecosm/apecosm_orca1/processed_pacific/obs_equatorial_mean.nc')
