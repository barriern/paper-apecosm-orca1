# +
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
import sys
sys.path.append('../../nino')
from extract_nino import read_index
import apecosm.ts as ts
import scipy.signal as sig
plt.rcParams['image.cmap'] = 'RdBu_r'

letters = list(string.ascii_lowercase)
letters = letters[1:]

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

gridparams = {'crs': ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}

mesh = xr.open_dataset('../../data/mesh_mask_eORCA1_v2.2.nc')
mesh = mesh.isel(t=0, z=0)
lon = mesh['glamt'].values
lat = mesh['gphit'].values
lonf = mesh['glamf'].values
latf = mesh['gphif'].values
tmask = mesh['tmask'].values

fig = plt.figure(figsize=(8, 16), facecolor='white')
plt.subplots_adjust(top=0.95)
axes_class = (GeoAxes, dict(map_projection=proj))

left = 0.06
width = 1 - 2 * left
bottom = 0.1
height = 0.45

axes = (left, bottom, width, height)
axgr = AxesGrid(fig, axes,  axes_class=axes_class, nrows_ncols=(2, 2), axes_pad=(0.7, 0.75), label_mode='', cbar_mode='each', cbar_size=0.1, cbar_pad=0.3, cbar_location="bottom")

axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

################################################################### Plotting covariances OBS

hadley = xr.open_dataset('data/cov_hadsst_oni_tpi.nc')
lonhad = hadley['lon'].values
lathad = hadley['lat'].values
hadoni = hadley['covoni'].to_masked_array().T
hadoni = np.ma.masked_where(hadoni == 0, hadoni)

model = xr.open_dataset('data/cov_modssh_oni_tpi.nc')
modsshoni = model['covoni'].to_masked_array().T * 100
modsshoni = np.ma.masked_where(modsshoni == 0, modsshoni)

model = xr.open_dataset('data/cov_modsst_oni_tpi.nc')
modoni = model['covoni'].to_masked_array().T
modoni = np.ma.masked_where(modoni == 0, modoni)

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
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])
ax.set_ylim(-40, 40)
ax.set_xlim(-60, 130)
ax.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

ax.set_title('Hadley SST / ONI')
cbax = axgr.cbar_axes[iiii]
cb = cbax.colorbar(cs)
try:
    cb.set_label('Covariance [C]')
except:
    cb.set_label_text('Covariance [C]')
    
################################################################### Plotting covariances SST MOD
    
iiii = 2
ax = axgr[iiii]

cs = ax.pcolormesh(lonf, latf, modoni[1:, 1:], transform=proj2)
cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax.add_feature(cfeature.COASTLINE, zorder=1001)
ax.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

gl = ax.gridlines(**gridparams)
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
xticks = np.array([150, 180, -180, -150, -120, -90, -60])
gl.xlocator = mticker.FixedLocator(xticks)
ax.set_ylim(-40, 40)
ax.set_xlim(-60, 130)

cbax = axgr.cbar_axes[iiii]
cb = cbax.colorbar(cs)
try:
    cb.set_label('Covariance [C]')
except:
    cb.set_label_text('Covariance [C]')

ax.set_title('Model SST / ONI')

################################################################### Plotting covariances SSH MOD
    
iiii = 1
ax = axgr[iiii]

cs = ax.pcolormesh(lonf, latf, modsshoni[1:, 1:], transform=proj2)
#cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax.add_feature(cfeature.COASTLINE, zorder=1001)
ax.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

gl = ax.gridlines(**gridparams)
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
xticks = np.array([150, 180, -180, -150, -120, -90, -60])
gl.xlocator = mticker.FixedLocator(xticks)
ax.set_ylim(-40, 40)
ax.set_xlim(-60, 130)
cs.set_clim(-8, 8)

cbax = axgr.cbar_axes[iiii]
cb = cbax.colorbar(cs)
try:
    cb.set_label('Covariance [cm]')
except:
    cb.set_label_text('Covariance [cm]')

ax.set_title('Model SSH / ONI')


################################################################### Plotting time-series ONI index

dnino, nino = read_index('../../data/index/oni.data')
iok = np.nonzero((dnino >= 195801) & (dnino <= 201812))[0]
dnino = dnino[iok]
nino = nino[iok]
time = np.arange(len(dnino))
ynino = dnino // 100
mnino = dnino - 100 * ynino

labels = np.array(['%.4d-%.2d' %(y, m) for y,m in zip(ynino, mnino)])

data = xr.open_dataset("data/simulated_enso_index.nc")
years = data['time_counter.year'].values 
months = data['time_counter.month'].values
enso = data['enso'].rolling(time_counter=3, center=True).mean()
#enso = data['enso'].values
nmod = len(time)

istart = np.nonzero(dnino == 195802)[0][0]
iend = np.nonzero(dnino == 201811)[0][0] + 1

test = np.corrcoef(enso[1:-1], nino[istart:iend])[0, 1]
print('Correlation ONI', test)

istart = np.nonzero(dnino == 195801)[0][0]
iend = np.nonzero(dnino == 201812)[0][0] + 1
test = np.corrcoef(enso, nino[istart:iend])

xticks = np.arange(2*12, len(time), 5 * 12)

left = 0.07
width = 0.85
bottom = 0.55
height = 0.1

axes = (left, bottom, width, height)

alpha = 0.7

ax = plt.axes(axes)
l1 = plt.fill_between(time, 0, nino, where=(nino>0), color='firebrick', interpolate=True)
l2 = plt.fill_between(time, 0, nino, where=(nino<0), color='steelblue', interpolate=True)
l3 = plt.plot(time, enso, 'k', label='Sim.', alpha=alpha)
plt.legend([l1, l3[0]], ['Obs.', 'Model'], loc=0, fontsize=8, ncol=2)
ax.set_title('ONI')

#plt.legend(loc=0)
ax.set_xticks(time[xticks])
ax.set_xticklabels(labels[xticks], rotation=45, ha='right')
ax.grid(True)
ax.set_xlim(time.min(), time.max())
ax.set_ylim(-4, 4)
ax.text(time[-1] - 50, -3, 'a' + ")", ha='center', va='center', bbox=dicttext)

plt.savefig('fig1', bbox_inches='tight')
# -


