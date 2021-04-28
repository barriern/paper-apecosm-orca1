import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import string

letters = list(string.ascii_lowercase)

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

gridparams = {'crs': ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}

mesh = xr.open_dataset('../data/mesh_mask_eORCA1_v2.2.nc')
mesh = mesh.isel(t=0, z=0)
lon = mesh['glamt'].values
lat = mesh['gphit'].values
lonf = mesh['glamf'].values
latf = mesh['gphif'].values
tmask = mesh['tmask'].values

hadley = xr.open_dataset('cov_hadsst_oni_tpi.nc')
lonhad = hadley['lon'].values
lathad = hadley['lat'].values
hadoni = hadley['covoni'].to_masked_array().T
hadoni = np.ma.masked_where(hadoni == 0, hadoni)
hadtpi = hadley['covtpi'].to_masked_array().T
hadtpi = np.ma.masked_where(hadtpi == 0, hadtpi)

model = xr.open_dataset('cov_modsst_oni_tpi.nc')
modoni = model['covoni'].to_masked_array().T
modoni = np.ma.masked_where(modoni == 0, modoni)
modtpi = model['covtpi'].to_masked_array().T
modtpi = np.ma.masked_where(modtpi == 0, modtpi)


fig = plt.figure(figsize=(12, 8), dpi=200)
axes_class = (GeoAxes, dict(map_projection=proj))
axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(2, 2), axes_pad=(0.7, 0.7), label_mode='', cbar_mode='each', cbar_size=0.1, cbar_pad=0.3, cbar_location="bottom")

axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

lontext = 120
lattext = 30
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)

iiii = 0
ax = axgr[iiii]
ccc = 1.1
step = 0.1
levels = np.arange(-ccc, ccc + step, step)

cs = ax.pcolormesh(lonhad, lathad, hadoni, transform=proj2, shading='auto')
cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax.add_feature(cfeature.COASTLINE, zorder=1001)

gl = ax.gridlines(**gridparams)
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])
ax.set_ylim(-40, 40)
ax.set_xlim(-60, 130)
ax.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

ax.set_title('Hadley SST')
cbax = axgr.cbar_axes[iiii]
cb = cbax.colorbar(cs)
cb.set_label('SST covariance (C)')
iiii = 1 
ax = axgr[iiii]

cs = ax.pcolormesh(lonf, latf, modoni[1:, 1:], transform=proj2)
cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax.add_feature(cfeature.COASTLINE, zorder=1001)
ax.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

gl = ax.gridlines(**gridparams)
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
xticks = np.array([150, 180, -180, -150, -120, -90, -60])
gl.xlocator = mticker.FixedLocator(xticks)
ax.set_ylim(-40, 40)
ax.set_xlim(-60, 130)

cbax = axgr.cbar_axes[iiii]
cb = cbax.colorbar(cs)
cb.set_label('SST covariance (C)')

ax.set_title('Model SST')

iiii = 2
ax = axgr[iiii]
ccc = 0.1
step = 0.1
levels = np.arange(-ccc, ccc + step, step)

cs = ax.pcolormesh(lonhad, lathad, hadtpi, transform=proj2, shading='auto')
cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax.add_feature(cfeature.COASTLINE, zorder=1001)
ax.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

gl = ax.gridlines(**gridparams)
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])
ax.set_ylim(-40, 40)
ax.set_xlim(-60, 130)

ax.set_title('Hadley SST')
cbax = axgr.cbar_axes[iiii]
cb = cbax.colorbar(cs)
cb.set_label('SST covariance (C)')

iiii = 3
ax = axgr[iiii]

cs = ax.pcolormesh(lonf, latf, modtpi[1:, 1:], transform=proj2)
cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax.add_feature(cfeature.COASTLINE, zorder=1001)

gl = ax.gridlines(**gridparams)
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
xticks = np.array([150, 180, -180, -150, -120, -90, -60])
gl.xlocator = mticker.FixedLocator(xticks)
ax.set_ylim(-40, 40)
ax.set_xlim(-60, 130)
ax.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

cbax = axgr.cbar_axes[iiii]
cb = cbax.colorbar(cs)
cb.set_label('SST covariance (C)')

ax.set_title('Model SST')

plt.savefig('covariance_maps_hadley_model.png', bbox_inches='tight')
