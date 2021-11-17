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
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc').isel(gpred=0)
const = const.rename({'wpred': 'w'})
const = const.assign(w=const['length'].values * 100)
const

ilat = slice(None, -3)
# -

mesh = xr.open_dataset('../new-97/data/pacific_mesh_mask.nc').isel(z=0, y=ilat)
lonf = mesh['glamf'].values
latf = mesh['gphif'].values

data = xr.open_dataset('../new-97/data/pacific_nino97_OOPE.nc').isel(y=ilat)
data = data.assign(w=const['length'] * 100)
data

tot = data['OOPE'] * const['weight_step']
tot = tot.sel(w=slice(25, 75))
tot = tot.sum(dim='w')
tot


def extract_ond(year):
    
    datestart = '%.4d-10-01' %year
    datestart

    dateend = '%.4d-12-31' %year
    dateend

    catch_96 = tot.sel(time=slice(datestart, dateend))
    output = catch_96.mean(dim='time')
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

axgr = ImageGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(3, 1), axes_pad=(0.0, 0.6), label_mode='',
                cbar_mode='each', cbar_size=0.1, cbar_pad=0.3, cbar_location="right", share_all=True)

diff = catch_97 - catch_96
tp_catch_97 = catch_97.where(catch_97 != 0)
tp_catch_96 = catch_96.where(catch_96 != 0)
diff = diff.where(diff != 0)

gridparams = {'crs':ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 
              'color':'k', 'linestyle':'--'}

cmax = 500
cmin = 100
ccc = 500
#cmin = 0
#cmax = 3

ax = axgr[0]
cs = ax.pcolormesh(lonf, latf, tp_catch_96.values[1:, 1:], transform=projin, cmap=plt.cm.jet)
cb = plt.colorbar(cs, axgr.cbar_axes[0])
cb.set_label('Catch (MT)')
ax.add_feature(cfeature.LAND, zorder=2)
ax.add_feature(cfeature.COASTLINE, zorder=3)
ax.set_extent([120.5, 299.5, -39.5, 39.5], crs=projin)
cs.set_clim(cmin, cmax)
ax.set_title('OND 96')
gl = ax.gridlines(**gridparams)
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])


ax = axgr[1]
cs = ax.pcolormesh(lonf, latf, tp_catch_97.values[1:, 1:], transform=projin, cmap=plt.cm.jet)
cb = plt.colorbar(cs, axgr.cbar_axes[1])
cb.set_label('Catch (MT)')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
cs.set_clim(cmin, cmax)
ax.set_title('OND 97')
gl = ax.gridlines(**gridparams)
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])

ax = axgr[2]
cs = ax.pcolormesh(lonf, latf, diff.values[1:, 1:], transform=projin, shading='auto')
cb = plt.colorbar(cs, axgr.cbar_axes[2])
cb.set_label('Catch (MT)')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
cs.set_clim(-ccc, ccc)
ax.set_title('OND 97 - OND 96')
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
plt.savefig('plot_ond_sardara_catches.png', bbox_inches='tight')
gl = ax.gridlines(**gridparams)
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])

plt.savefig('plot_ond_oope_diff.png', bbox_inches='tight')
# -


