import numpy as np
from datetime import date
import os.path
import numpy as np
from eofs.standard import Eof
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import stats
import xarray as xr
from glob import glob
import apecosm.ts as ts

dirout = './'

# Load the mesh mask
mesh = xr.open_dataset("../../data/mesh_mask_eORCA1_v2.2.nc")
tmask = np.squeeze(mesh['tmask'].values[0, 0])
e1t = mesh['e1t'].values
e2t = mesh['e2t'].values
surf = e1t * e2t
surf = np.squeeze(surf)
nlat, nlon = surf.shape
lon = mesh['glamt'].values
lat = mesh['gphit'].values
lon = np.squeeze(lon)
lat = np.squeeze(lat)
lonf = mesh['glamf'].values[0]
latf = mesh['gphif'].values[0]

latmin = -5
latmax = 5
lonmax = -120
lonmin = -170

test = (lat<=latmax) & (lat>=latmin)
test = test & (lon<=lonmax) & (lon>=lonmin)
test = test & (tmask == 1)

ilat, ilon = np.nonzero(test == True)

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

output = tmask.copy()
output[ilat, ilon] = 2

plt.figure()
ax = plt.gca(projection=proj)
cs = ax.pcolormesh(lonf, latf, output[1:, 1:], transform=proj2)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
plt.colorbar(cs)
plt.savefig('boxes.png')


data = xr.open_mfdataset('data/nemo/*nc')
years = data['time_counter.year'].values
months = data['time_counter.month'].values
date = years * 100 + months

datestart = 197101
dateend = 200012

iclim = np.nonzero((date >= datestart) & (date <= dateend))[0]

sst = np.squeeze(data['thetao'].values)
ntime, nlat, nlon = sst.shape

clim, anom = ts.get_monthly_clim(sst[iclim, :, :])

nyears = ntime // 12
index = np.arange(12)

anom = np.zeros(sst.shape)
for y in range(nyears):
    anom[index] = sst[index] - clim
    index += 12

surf = surf[np.newaxis, :, :] 

numer = np.sum(surf[:, ilat, ilon] * anom[:, ilat, ilon], axis=-1)
denom = np.sum(surf[:, ilat, ilon], axis=-1)
output = numer / denom

timeout = data['time_counter.year'].values * 100 + data['time_counter.month'].values 

dsout = xr.Dataset({'enso':(['time'], output)})
dsout['time'] = (['time'], timeout)
dsout.attrs['script'] = os.path.realpath(__file__)
dsout.attrs['description'] = 'EOF for densities by class'
dsout.to_netcdf('%s/simulated_enso_index.nc' %(dirout))
