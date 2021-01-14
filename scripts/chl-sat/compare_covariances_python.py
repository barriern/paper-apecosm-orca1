import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

gridparams = {'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}

#gl = ax.gridlines(crs=proj, draw_labels=True,
#                  linewidth=2, color='gray', alpha=0.5, linestyle='--')

plt.figure()
plt.subplots_adjust(top=0.95, bottom=0.2, hspace=0.25)
ccc = 0.05
space = 0.01
levels = np.arange(-0.1, 0.1 + space, space)

# processing obs
data = xr.open_dataset("interp/covariance_satellite_data.nc")
cov = data['cov'].to_masked_array()
cov = np.ma.masked_where(cov == 0, cov)
lon = data['lon'].values
lat = data['lat'].values

ax = plt.subplot(2, 1, 1, projection=proj)
ax.set_ylim(-40, 40)
ax.set_xlim(-60, 130)

cs = plt.pcolormesh(lon, lat, cov, transform=ccrs.PlateCarree())
cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax.add_feature(cfeature.COASTLINE, zorder=1001)

gl = ax.gridlines(**gridparams)
gl.xlabels_top = False
gl.ylabels_right = False
#gl.xlines = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -150, -120, -90, -60])

plt.title('Obs.')

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

ax2 = plt.subplot(2, 1, 2, projection=proj, sharex=ax, sharey=ax)

cs = plt.pcolormesh(lon, lat, cov, transform=proj2)
cs.set_clim(-ccc, ccc)
ax2.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax2.add_feature(cfeature.COASTLINE, zorder=1001)

plt.title('Model')

xmin = 0.2
cax = plt.axes([xmin, 0.1, 1-2*xmin, 0.03])
cb = plt.colorbar(cs, cax, orientation='horizontal')
cb.set_label("Chl. cov (mg/m3)")

gl = ax2.gridlines(**gridparams)
gl.xlabels_top = False
gl.ylabels_right = False
#gl.ylabels_left = False
#gl.xlines = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -150, -120, -90, -60])

plt.savefig('compare_covariance_chl.png', bbox_inches='tight')
