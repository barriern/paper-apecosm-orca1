# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import cftime
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

species = 'SKJ'
gear = 'PS'
res = 5

data = xr.open_dataset('data/regridded_catch_gear_%s_species_%s_%dx%d.nc' %(gear, species, res, res), decode_times=False)
data
# -

# ## Date corrections
#
# Dates are converted into datetime object

dates = data['time'].values
dates[:5]

years = dates // 100
years[:5]

months = dates - 100 * years
months[:5]

time = [datetime(y, m, 1) for y, m in zip(years, months)]
time[:5]

data['time'] = time

# ## Masking the missing data

data = data.where(data['catch'] != 0)
data

# ## Computation of the climatology

dataclim = data.sel(time=slice('1970-01-01', '2000-12-31'))
dataclim

values = dataclim['catch'].to_masked_array()

# +
ntime, nlat, nlon = values.shape
nyears = ntime // 12
index = np.arange(12)

clim = np.ma.zeros((12, nlat, nlon))
count = np.zeros((12, nlat, nlon))

for i in range(nyears):
    
    temp = values[index, :, :]
    
    mask = np.ma.getmaskarray(values[index, :, :])  # 1 if mask, 0 else

    clim[mask == 0] += temp[mask == 0]
    
    count += (1 - mask).astype(int)
    
    index += 12
    
clim /= nyears
count.shape
# -

count = (np.sum(count, axis=0) / ntime) * 100
count.shape

clim = np.mean(clim, axis=0)
clim.shape

clim = np.ma.masked_where(count == 0, clim)

# ## Normalisation by the surface

lat = data['lat'].values
lat

Rt = 6371.009
surface = Rt * np.deg2rad(res) * Rt * np.deg2rad(res) * np.cos(np.deg2rad(lat))

clim = clim / surface[:, np.newaxis]

# +
fig = plt.figure(figsize=(12, 12))
plt.rcParams['contour.negative_linestyle'] = 'solid'
plt.rcParams['font.size'] = 15
thres = 30

lon = data['lon'].values
lat = data['lat'].values
projin = ccrs.PlateCarree()
projout = ccrs.PlateCarree(central_longitude=180)
ax = plt.axes(projection=projout)
ax.coastlines()
ax.add_feature(cfeature.LAND)
cs = ax.pcolormesh(lon, lat, np.log10(clim, where=clim.mask == False), shading='auto', transform=projin)
cb = plt.colorbar(cs, orientation='horizontal', pad=0.02)
ax.set_title('Mean catch, species=%s, gear=%s' %(species, gear))
cb.set_label('MT/km2')
iok = np.nonzero(count >= thres)
plt.plot(lon[iok[1]], lat[iok[0]], marker='.', linestyle='none', transform=projin)
plt.savefig('mean_catch_gear_%s_species_%s_%dx%d.png' %(gear, species, res, res))
# -


