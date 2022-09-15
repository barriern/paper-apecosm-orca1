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

mesh = mesh.where(abs(lat) == 0)
mesh['tmask'].isel(z=0).plot()

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
# -

boolean = (lon0.values < 270) & (lon0.values >= 150)
iok = np.nonzero(boolean)[0]
iok

mesh['tmask'].isel(z=0, x=iok).plot()

clim = xr.open_dataset('data/mean_thetao.nc')
iok = np.nonzero(lon0.values == 150.5)[0][0]
iok

# +
fig = plt.figure(figsize=(16, 12))

off = 1

axgr = []
for i in range(6):
    axgr.append(plt.subplot(3, 2, i + 1))

plt.rcParams['lines.linewidth'] = 0.5
plt.rcParams['font.size'] = 17
plt.rcParams['image.cmap'] = 'RdBu_r'

lontext = 260
ztext = -210
lontext2 = 160

dictext2 = dict(ha='left', va='center', bbox=dicttext, zorder=20)

# thetao (0)
# uo (1)
# plk (2)

# plottiing thetao
cpt = 0
ax = axgr[cpt]
clim = xr.open_dataset('data/mean_thetao.nc').isel(x=iok)
clim = clim['__xarray_dataarray_variable__'].to_masked_array()

anom = xr.open_dataset('data/mean_thetao_anomalies_ond_97.nc').isel(x=iok)
anom = anom['__xarray_dataarray_variable__'].values

cl = ax.plot(clim, -depth[:], )
cs = ax.plot(anom + off * clim, -depth[:])
plt.title('temp')

# plotiing uo
cpt = 2
ax = axgr[cpt]
clim = xr.open_dataset('data/mean_uo.nc').isel(x=iok)
clim = clim['__xarray_dataarray_variable__'].to_masked_array()

anom = xr.open_dataset('data/mean_uo_anomalies_ond_97.nc').isel(x=iok)
anom = anom['__xarray_dataarray_variable__'].to_masked_array()

step = 0.1
cl = ax.plot(clim[:-1], -depth[:-1], )
cs = ax.plot((anom + off * clim)[:-1], -depth[:-1],)
ax.set_ylabel('Depth (m)')
plt.title('U')

# plotiing PLK
cpt = 4
ax = axgr[cpt]
clim = xr.open_dataset('data/mean_PLK.nc').isel(x=iok)
clim = clim['__xarray_dataarray_variable__'].to_masked_array()

anom = xr.open_dataset('data/mean_PLK_anomalies_ond_97.nc').isel(x=iok)
anom = anom['__xarray_dataarray_variable__'].to_masked_array()
step = 0.05
cl = ax.plot( clim[:-1], -depth[:-1],)
cs = ax.plot((anom + off * clim)[:-1], -depth[:-1], )
ax.set_ylabel('Depth (m)')
plt.title('plk')

# ### plotting forage
clim = xr.open_dataset('data/mean_forage.nc').isel(x=iok)
clim = clim['__xarray_dataarray_variable__'].to_masked_array()

anom = xr.open_dataset('data/mean_forage_anomalies_ond_97.nc').isel(x=iok)
anom = anom['__xarray_dataarray_variable__'].to_masked_array()

print(clim.shape)

cpt = 1
ax = axgr[cpt]
step = 0.05
#cl = ax.plot(clim[:, 0].T[:-1],  -depth[:-1], )
cs = ax.plot((anom + off * clim)[:, 0].T[:-1], -depth[:-1],)

cpt = 3
ax = axgr[cpt]
step = 0.05
#cl = ax.plot(clim[:, 1].T[:-1], -depth[:-1])
cs = ax.plot((anom + off * clim)[:, 1].T[:-1], -depth[:-1])

cpt = 5
ax = axgr[cpt]
step = 0.05
#cl = ax.plot(clim[:, 2].T[:-1], -depth[:-1])
cs = ax.plot((anom + off * clim)[:, 2].T[:-1], -depth[:-1])
#labels = ['150', '180', '-150', '-120', '-90']
#xticks = np.array([float(l) for l in labels])
#xticks[xticks < 0] += 360
#ax.set_xticks(xticks)

plt.savefig('forage_mean_ond97.png', bbox_inches='tight', facecolor='white')
# -

