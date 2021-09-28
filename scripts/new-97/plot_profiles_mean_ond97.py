# +
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
import numpy as np

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


# +
fig = plt.figure(figsize=(14, 8))
axgr = ImageGrid(fig, 111, nrows_ncols=(3, 1), ngrids=None, direction='row', axes_pad=(0.8, 0.2), share_all=False, aspect=True, label_mode='L', cbar_mode='each', cbar_location='right', cbar_pad='5%', cbar_size='5%', cbar_set_cax=True, axes_class=None)
cbar_axes = axgr.cbar_axes

plt.rcParams['lines.linewidth'] = 0.5

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
cb.set_label('Temp. (C)')
cs.set_clim(-8, 8)

# plotiing uo
cpt = 1
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
cb.set_label('Zonal Vel. (m/s)')

cs.set_clim(-0.6, 0.6)

# plotiing PLK
cpt = 2
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
cb.set_label('LTL conc. (mmol/m3)')

#cs.set_clim(-0.6, 0.6)
# -


