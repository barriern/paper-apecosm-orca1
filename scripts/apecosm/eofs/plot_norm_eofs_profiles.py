import matplotlib
import xarray as xr
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
import matplotlib

signs = -np.ones((100))
index = list(range(84, 93 + 1)) + [99]
index = np.array(index)
signs[index] *= -1

matplotlib.rcParams['image.cmap'] = 'RdBu_r'

const = xr.open_dataset('../../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
length = const['length'].values * 100
wstep = const['weight_step'].values

ilength = np.nonzero(length <= 100)[0]
length = length[ilength]

mesh = xr.open_dataset("../../data/mesh_mask_eORCA1_v2.2.nc")
mesh = mesh.isel(x=slice(58, 228, None), y=slice(140, 233, None))
lon = mesh['glamt'].values[0]
lat = mesh['gphit'].values[0]

data = xr.open_dataset('data/eof_full_density_20.nc')

eofmap0 = data['eofmap0'].to_masked_array()  # bins, eofs, lat, lon
eofmap0 = (eofmap0.T * signs).T
eofmap0 = eofmap0[ilength, 0, :, :]  # size, lat, lon

lon = np.mean(lon, axis=0)
lat = np.mean(lat, axis=1)

lon0 = -150
ilon = np.argmin((lon - lon0)**2)
print(ilon, lon[ilon])

ilat = np.nonzero(np.abs(lat) <= 20)[0]
lat = lat[ilat]
eofmap0 = eofmap0[:, ilat, :]  

cc0 = 0.04
step = 0.01
levels = np.arange(-cc0, cc0 + step, step)

plt.figure()
plt.subplots_adjust(hspace=0.4)
ax = plt.subplot(2, 2, 1)
tp = eofmap0[:, :, ilon].T
cc0 = np.percentile(np.abs(tp[tp.mask == False]), 99.9)
cs = ax.pcolormesh(length, lat, eofmap0[:, :, ilon].T)
cl = ax.contour(length, lat, eofmap0[:, :, ilon].T, levels=levels, colors='k', linewidths=0.5)
cl2 = ax.contour(length, lat, eofmap0[:, :, ilon].T, levels=0, colors='k', linewidths=1)
cs.set_clim(-cc0, cc0)
cb = plt.colorbar(cs)
cb.add_lines(cl)
ax.set_xscale('log')
#ax.set_xlabel('Length (cm)')
ax.set_ylabel('Latitude')
plt.title('150W')
ax.grid(True, linestyle='--', color='gray', linewidth=0.5)

lon0 = 160
ilon = np.argmin((lon - lon0)**2)
print(ilon, lon[ilon])

ax = plt.subplot(2, 2, 2)
tp = eofmap0[:, :, ilon].T
cc0 = np.percentile(np.abs(tp[tp.mask == False]), 99.9)
cs = ax.pcolormesh(length, lat, eofmap0[:, :, ilon].T)
cl = ax.contour(length, lat, eofmap0[:, :, ilon].T, levels=levels, colors='k', linewidths=0.5)
cl2 = ax.contour(length, lat, eofmap0[:, :, ilon].T, levels=0, colors='k', linewidths=1)
cs.set_clim(-cc0, cc0)
cb = plt.colorbar(cs)
cb.add_lines(cl)
ax.set_xscale('log')
plt.title('160E')
ax.set_xlabel('Length (cm)')
ax.set_ylabel('Latitude')
ax.grid(True, linestyle='--', color='gray', linewidth=0.5)

lon0 = -120
ilon = np.argmin((lon - lon0)**2)
print(ilon, lon[ilon])

ax = plt.subplot(2, 2, 3)
tp = eofmap0[:, :, ilon].T
cc0 = np.percentile(np.abs(tp[tp.mask == False]), 99.9)
print(cc0)
cs = ax.pcolormesh(length, lat, eofmap0[:, :, ilon].T)
cl = ax.contour(length, lat, eofmap0[:, :, ilon].T, levels=levels, colors='k', linewidths=0.5)
cl2 = ax.contour(length, lat, eofmap0[:, :, ilon].T, levels=0, colors='k', linewidths=1)
cs.set_clim(-cc0, cc0)
cb = plt.colorbar(cs)
cb.add_lines(cl)
ax.set_xscale('log')
plt.title('120W')
ax.set_xlabel('Length (cm)')
ax.set_ylabel('Latitude')
ax.grid(True, linestyle='--', color='gray', linewidth=0.5)

plt.savefig('profiles_eofs.pdf')

