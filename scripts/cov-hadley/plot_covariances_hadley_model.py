import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes

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

hadley = xr.open_dataset('data//covariance_yearly_enso_surface_sst.nc')
lonhad = hadley['longitude'].values
lathad = hadley['latitude'].values
hadley = hadley['covariance'].to_masked_array()
hadley = np.ma.masked_where(hadley == 0, hadley)

model = xr.open_dataset('data//covariance_yearly_enso_surface_thetao.nc')
model = model['covariance'].to_masked_array()
model = np.ma.masked_where(model == 0, model)

fig = plt.figure(figsize=(12, 8), dpi=200)
axes_class = (GeoAxes, dict(map_projection=proj))
axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(2, 1), axes_pad=(1, 0.5), label_mode='', cbar_mode='single', cbar_size=0.1, cbar_pad=0., cbar_location="bottom")
ccc = 1.5

axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

ax = axgr[0]
ccc = 1.5
step = 0.1
levels = np.arange(-ccc, ccc + step, step)

cs = ax.pcolormesh(lonhad, lathad, hadley, transform=proj2, shading='auto')
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

ax.set_title('Hadley SST')

ax = axgr[1]

cs = ax.pcolormesh(lonf, latf, model[1:, 1:], transform=proj2)
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

cbax = axgr.cbar_axes[0]
cb = cbax.colorbar(cs)
cb.set_label('SST covariance (C)')

ax.set_title('Model SST')

plt.savefig('covariance_maps_hadley_model.png', bbox_inches='tight')
