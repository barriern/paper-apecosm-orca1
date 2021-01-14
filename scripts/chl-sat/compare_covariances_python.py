import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

gridparams = {'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}

#gl = ax.gridlines(crs=proj, draw_labels=True,
#                  linewidth=2, color='gray', alpha=0.5, linestyle='--')

plt.figure(figsize=(12, 8))
fig = plt.figure(figsize=(12, 8))
axes_class = (GeoAxes, dict(map_projection=proj))
axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(2, 1), axes_pad=(1, 0.5), label_mode='', cbar_mode='single', cbar_size=0.1, cbar_pad=0., cbar_location="bottom", share_all=True)

ccc = 0.05
space = 0.01
levels = np.arange(-0.1, 0.1 + space, space)

# processing obs
data = xr.open_dataset("interp/covariance_satellite_data.nc")
cov = data['cov'].to_masked_array()
cov = np.ma.masked_where(cov == 0, cov)
lon = data['lon'].values
lat = data['lat'].values

ax = axgr[0]
ax.set_ylim(-40, 40)
ax.set_xlim(-60, 130)

cs = ax.pcolormesh(lon, lat, cov, transform=ccrs.PlateCarree())
cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax.add_feature(cfeature.COASTLINE, zorder=1001)

gl = ax.gridlines(**gridparams)
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])

ax.set_title('Obs.')

#######

mesh = xr.open_dataset('/home/barrier/Work/apecosm/ORCA1/DATA_APE_ORCA1/pisces/data/corrected_mesh_mask_eORCA1_v2.2.nc')
mesh = mesh.isel(t=0)
lon = mesh['glamt'].values
lat = mesh['gphit'].values
tmask = mesh['tmask'].values[0]

data = xr.open_dataset("model/covariance_model_data.nc")
cov = data['cov'].values
cov = np.ma.masked_where(tmask == 0, cov)
print(cov.min(), cov.max())

ax2 = axgr[1]

cs = ax2.pcolormesh(lon, lat, cov, transform=proj2)
cs.set_clim(-ccc, ccc)
ax2.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax2.add_feature(cfeature.COASTLINE, zorder=1001)

ax2.set_title('Model')

xmin = 0.2
#cax = plt.axes([xmin, 0.1, 1-2*xmin, 0.03])
cb = axgr.cbar_axes[0].colorbar(cs)
cb.set_label_text("Chl. cov (mg/m3)")

gl = ax2.gridlines(**gridparams)
gl.xlabels_top = False
gl.ylabels_right = False
#gl.ylabels_left = False
#gl.xlines = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])

plt.savefig('compare_covariance_chl.png', bbox_inches='tight')
