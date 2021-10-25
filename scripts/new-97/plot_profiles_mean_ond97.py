# +
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
import numpy as np
from cartopy.mpl.ticker import (LatitudeFormatter, LongitudeFormatter,
                                LatitudeLocator, LongitudeLocator)
from cartopy.mpl.geoaxes import GeoAxes
import cartopy.crs as ccrs

import string
letters = list(string.ascii_lowercase)
# -

mesh = xr.open_dataset('data/pacific_mesh_mask.nc')

lat = mesh['gphit']
lat

mesh = mesh.where(abs(lat) <= 5)
mesh

lon0 = mesh['glamt'].mean(dim='y').isel(z=0)
lon0

lon0 = (lon0 + 360) % 360

depth = mesh['gdept_1d'].mean(dim='y').isel(x=0)
depth

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)


# Now, we define formatters (did not find a better way to do it...)

# +
# central 0
formatter0 = LongitudeFormatter(dateline_direction_label=True)
p = ccrs.PlateCarree()
test_ticks = [-180, -120, -60, 0, 60, 120, 180]
result = [formatter0(tick) for tick in test_ticks]
expected = ['180\u00B0W', '120\u00B0W', '60\u00B0W', '0\u00B0',
            '60\u00B0E', '120\u00B0E', '180\u00B0E']
print(result)
print(expected)

# central 180
formatter180 = LongitudeFormatter(zero_direction_label=True)
p = ccrs.PlateCarree(central_longitude=180)
temp = plt.axes(projection=p, frame_on=False)
formatter180.axis = temp
del(temp)
test_ticks = [-180, -120, -60, 0, 60, 120, 180]
result = [formatter180(tick) for tick in test_ticks]
expected = ['0\u00B0E', '60\u00B0E', '120\u00B0E', '180\u00B0',
            '120\u00B0W', '60\u00B0W', '0\u00B0W']
print(result)
print(expected)

# +
fig = plt.figure(figsize=(16, 12))
axgr = ImageGrid(fig, 111, nrows_ncols=(3, 2), ngrids=None, direction='row', axes_pad=(1, 0.2), share_all=True, aspect=False, label_mode='L', cbar_mode='each', cbar_location='right', cbar_pad='5%', cbar_size='5%', cbar_set_cax=True, axes_class=None)
cbar_axes = axgr.cbar_axes

plt.rcParams['lines.linewidth'] = 0.5
plt.rcParams['font.size'] = 17

# thetao (0)
# uo (1)
# plk (2)

# plottiing thetao
cpt = 0
ax = axgr[cpt]
clim = xr.open_dataset('data/mean_thetao.nc')
clim = clim['__xarray_dataarray_variable__'].to_masked_array()

anom = xr.open_dataset('data/mean_thetao_anomalies_ond_97.nc')
anom = anom['__xarray_dataarray_variable__'].to_masked_array()

cl = ax.contour(lon0, -depth[:-1], clim[:-1, :], 11, colors='k')
plt.clabel(cl)
cs = ax.pcolormesh(lon0, -depth[:-1], anom[:-1, :], shading='auto')
cb = plt.colorbar(cs, cbar_axes[cpt])
cb.set_label('[C]')
cs.set_clim(-8, 8)
ax.text(lontext2, ztext, 'T', ha='center', va='center', bbox=dicttext)
ax.text(lontext, ztext, 'a)', ha='center', va='center', bbox=dicttext)

lontext = 280
ztext = -210
lontext2 = 130

# plotiing uo
cpt = 2
ax = axgr[cpt]
clim = xr.open_dataset('data/mean_uo.nc')
clim = clim['__xarray_dataarray_variable__'].to_masked_array()

anom = xr.open_dataset('data/mean_uo_anomalies_ond_97.nc')
anom = anom['__xarray_dataarray_variable__'].to_masked_array()

step = 0.05
cl = ax.contour(lon0, -depth[:-1], clim[:-1, :], colors='k', levels = np.arange(-1, 1 + step, step))
plt.clabel(cl)
cs = ax.pcolormesh(lon0, -depth[:-1], anom[:-1, :], shading='auto')
cb = plt.colorbar(cs, cbar_axes[cpt])
cb.set_label('[m/s]')
ax.text(lontext2, ztext, 'U', ha='center', va='center', bbox=dicttext)
ax.text(lontext, ztext, 'b)', ha='center', va='center', bbox=dicttext)
ax.set_ylabel('Depth (m)')

cs.set_clim(-0.6, 0.6)

# plotiing PLK
cpt = 4
ax = axgr[cpt]
clim = xr.open_dataset('data/mean_PLK.nc')
clim = clim['__xarray_dataarray_variable__'].to_masked_array()

anom = xr.open_dataset('data/mean_PLK_anomalies_ond_97.nc')
anom = anom['__xarray_dataarray_variable__'].to_masked_array()
step = 0.05
cl = ax.contour(lon0, -depth[:-1], clim[:-1, :], 11, colors='k')
plt.clabel(cl)
cs = ax.pcolormesh(lon0, -depth[:-1], anom[:-1, :], shading='auto')
cb = plt.colorbar(cs, cbar_axes[cpt])
cb.set_label('[mmol/m3]')
print(anom.shape)
cs.set_clim(-2, 2)
ax.text(lontext2, ztext, 'LTL', ha='center', va='center', bbox=dicttext)
ax.text(lontext, ztext, 'c)', ha='center', va='center', bbox=dicttext)
ax.set_ylabel('Depth (m)')

### plotting forage
clim = xr.open_dataset('data/mean_forage.nc')
clim = clim['__xarray_dataarray_variable__'].to_masked_array()

anom = xr.open_dataset('data/mean_forage_anomalies_ond_97.nc')
anom = anom['__xarray_dataarray_variable__'].to_masked_array()

cpt = 1
ax = axgr[cpt]
step = 0.05
cl = ax.contour(lon0, -depth[:-1], clim[:, :, 0].T[:-1], 11, colors='k')
plt.clabel(cl)
cs = ax.pcolormesh(lon0, -depth[:-1], anom[:, :, 0].T[:-1], shading='auto')
cb = plt.colorbar(cs, cbar_axes[cpt])
cb.set_label('[J/m2]')
cs.set_clim(-50, 50)
ax.text(lontext2, ztext, '0-3cm', ha='center', va='center', bbox=dicttext)
ax.text(lontext, ztext, 'd)', ha='center', va='center', bbox=dicttext)
ax.set_ylabel('Depth (m)')

cpt = 3
ax = axgr[cpt]
step = 0.05
cl = ax.contour(lon0, -depth[:-1], clim[:, :, 1].T[:-1], 11, colors='k')
plt.clabel(cl)
cs = ax.pcolormesh(lon0, -depth[:-1], anom[:, :, 1].T[:-1], shading='auto')
cb = plt.colorbar(cs, cbar_axes[cpt])
cb.set_label('[J/m2]')
cs.set_clim(-20, 20)
ax.text(lontext2, ztext, '3-20cm', ha='center', va='center', bbox=dicttext)
ax.text(lontext, ztext, 'e)', ha='center', va='center', bbox=dicttext)
ax.set_ylabel('Depth (m)')

cpt = 5
ax = axgr[cpt]
step = 0.05
cl = ax.contour(lon0, -depth[:-1], clim[:, :, 2].T[:-1], 11, colors='k')
plt.clabel(cl)
cs = ax.pcolormesh(lon0, -depth[:-1], anom[:, :, 2].T[:-1], shading='auto')
cb = plt.colorbar(cs, cbar_axes[cpt])
cb.set_label('[J/m2]')
cs.set_clim(-6, 6)
ax.text(lontext2, ztext, '20-90cm', ha='center', va='center', bbox=dicttext)
ax.text(lontext, ztext, 'f)', ha='center', va='center', bbox=dicttext)
ax.set_ylabel('Depth (m)')
labels = ['150', '180', '-150', '-120', '-90', '-60']
xticks = np.array([float(l) for l in labels])
xticks[xticks < 0] += 360
ax.set_xticks(xticks)

ax.xaxis.set_major_formatter(formatter0)

plt.savefig('forage_mean_ond97.png', bbox_inches='tight', facecolor='white')
# -

