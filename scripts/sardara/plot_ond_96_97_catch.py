# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid
from cartopy.mpl.geoaxes import GeoAxes
import cartopy.feature as cfeature
import cartopy.crs as ccrs

skj = xr.open_dataset('data/regridded_catch_gear_PS_species_SKJ_1x1.nc')
skj
# -

yft = xr.open_dataset('data/regridded_catch_gear_PS_species_YFT_1x1.nc')
yft

tot = yft['catch'] + skj['catch']

dates = tot['time'].values


def extract_ond(year):
    
    datestart = year * 100 + 10
    datestart

    dateend = year * 100 + 12
    dateend

    iok = np.nonzero((dates >= datestart) & (dates <= dateend))[0]
    
    catch_96 = tot.isel(time=iok)
    output = catch_96.sum(dim='time')
    return output


catch_96 = extract_ond(1996)
catch_96

catch_97 = extract_ond(1997)
catch_97

# +
plt.rcParams['font.size'] = 15

projout = ccrs.PlateCarree(central_longitude=180)
projin = ccrs.PlateCarree()

fig = plt.figure(figsize=(8, 12), facecolor='white')
axes_class = (GeoAxes, dict(map_projection=projout))

axgr = ImageGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(3, 1), axes_pad=(0.7, 0.5), label_mode='',
                cbar_mode='each', cbar_size=0.1, cbar_pad=0.3, cbar_location="right", share_all=True)

lon = catch_96['lon'].values
lat = catch_96['lat'].values

diff = catch_97 - catch_96
tp_catch_97 = catch_97.where(catch_97 != 0)
tp_catch_96 = catch_96.where(catch_96 != 0)
diff = diff.where(diff != 0)

cmax = 1000
cmin = 100
ccc = 1000
#cmin = 0
#cmax = 3

ax = axgr[0]
cs = ax.pcolormesh(lon, lat, tp_catch_96.values, transform=projin, shading='auto', cmap=plt.cm.jet)
cb = plt.colorbar(cs, axgr.cbar_axes[0])
cb.set_label('Catch (MT)')
ax.add_feature(cfeature.LAND, zorder=2)
ax.add_feature(cfeature.COASTLINE, zorder=3)
ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=projin)
cs.set_clim(cmin, cmax)
ax.set_title('OND 96')

ax = axgr[1]
cs = ax.pcolormesh(lon, lat, tp_catch_97.values, transform=projin, shading='auto', cmap=plt.cm.jet)
cb = plt.colorbar(cs, axgr.cbar_axes[1])
cb.set_label('Catch (MT)')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
cs.set_clim(cmin, cmax)
ax.set_title('OND 97')

ax = axgr[2]
cs = ax.pcolormesh(lon, lat, diff.values, transform=projin, shading='auto')
cb = plt.colorbar(cs, axgr.cbar_axes[2])
cb.set_label('Catch (MT)')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
cs.set_clim(-ccc, ccc)
ax.set_title('OND 97 - OND 96')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
plt.savefig('plot_ond_sardara_catches.png', bbox_inches='tight')
# -


