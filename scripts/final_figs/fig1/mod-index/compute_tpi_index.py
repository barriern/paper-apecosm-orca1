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

''' Adapted from https://psl.noaa.gov/data/timeseries/IPOTPI/tpi.create.ncl'''

def extract_domain(latmin, latmax, lonmin, lonmax):
    test = ((lat >= latmin) & (lat <= latmax))
    test2 = (lon >= lonmin) | (lon <= lonmax)
    test = (test & test2) * tmask
    iok1 = np.nonzero(test)
    return iok1


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

output = tmask.copy()

latmin = 25
latmax= 45
lonmin = 140
lonmax = -145
iok1 = extract_domain(latmin, latmax, lonmin, lonmax)
output[iok1] = 100

latmin = -10
latmax= 10
lonmin = 170
lonmax = -90
iok2 = extract_domain(latmin, latmax, lonmin, lonmax)
output[iok2] = 200

latmin = -50
latmax= -15
lonmin = 150
lonmax = -160
iok3 = extract_domain(latmin, latmax, lonmin, lonmax)
output[iok3] = 300

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

plt.figure()
ax = plt.gca(projection=proj)
cs = ax.pcolormesh(lonf, latf, output[1:, 1:], transform=proj2)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
plt.colorbar(cs)
plt.savefig('boxes.png')

datestart = 197101
dateend = 200012

data = xr.open_mfdataset('data/nemo/*nc')
years = data['time_counter.year'].values
months = data['time_counter.month'].values
date = years * 100 + months

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

surf = surf[np.newaxis, :]

ilat = iok1[0]
ilon = iok1[1]
ts1 = np.sum(surf[:, ilat, ilon] * anom[:, ilat, ilon], axis=-1) / np.sum(surf[:, ilat, ilon], axis=-1)

ilat = iok2[0]
ilon = iok2[1]
ts2 = np.sum(surf[:, ilat, ilon] * anom[:, ilat, ilon], axis=-1) / np.sum(surf[:, ilat, ilon], axis=-1)

ilat = iok3[0]
ilon = iok3[1]
ts3 = np.sum(surf[:, ilat, ilon] * anom[:, ilat, ilon], axis=-1) / np.sum(surf[:, ilat, ilon], axis=-1)

tpi = ts2 - 0.5 * (ts1 + ts3)

dsout = xr.Dataset()
dsout['ts1'] = (['time'], ts1)
dsout['ts2'] = (['time'], ts2)
dsout['ts3'] = (['time'], ts3)
dsout['tpi'] = (['time'], tpi)
dsout['time'] = date
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('model_tpi_index.nc')
