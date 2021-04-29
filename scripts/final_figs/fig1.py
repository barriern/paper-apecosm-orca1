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
sys.path.append('../nino')
from extract_nino import read_index
import apecosm.ts as ts
import scipy.signal as sig
from envtoolkit.ts import Lanczos

nWgt      = 157
fca       = 1./156
fca = 1 / fca

lanc = Lanczos('lp', pca=fca, nwts=nWgt)


letters = list(string.ascii_lowercase)
letters = letters[2:]

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

fig = plt.figure(figsize=(8, 16))
plt.subplots_adjust(top=0.95)
axes_class = (GeoAxes, dict(map_projection=proj))

left = 0.06
width = 1 - 2 * left
bottom = 0.18 #0.5
height = 0.45

axes = (left, bottom, width, height)
axgr = AxesGrid(fig, axes,  axes_class=axes_class, nrows_ncols=(2, 2), axes_pad=(0.7, 0.75), label_mode='', cbar_mode='each', cbar_size=0.1, cbar_pad=0.3, cbar_location="bottom")

axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

#################################################################### Plotting covariances SST

hadley = xr.open_dataset('../cov-hadley/data/cov_hadsst_oni_tpi.nc')
lonhad = hadley['lon'].values
lathad = hadley['lat'].values
hadoni = hadley['covoni'].to_masked_array().T
hadoni = np.ma.masked_where(hadoni == 0, hadoni)
hadtpi = hadley['covtpi'].to_masked_array().T
hadtpi = np.ma.masked_where(hadtpi == 0, hadtpi)

model = xr.open_dataset('../cov-hadley/data/cov_modsst_oni_tpi.nc')
modoni = model['covoni'].to_masked_array().T
modoni = np.ma.masked_where(modoni == 0, modoni)
modtpi = model['covtpi'].to_masked_array().T
modtpi = np.ma.masked_where(modtpi == 0, modtpi)

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

ax.set_title('Hadley SST / ONI')
cbax = axgr.cbar_axes[iiii]
cb = cbax.colorbar(cs)
cb.set_label('Covariance [C]')
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
cb.set_label('Covariance [C]')

ax.set_title('Model SST / ONI')

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

ax.set_title('Hadley SST / TPI')
cbax = axgr.cbar_axes[iiii]
cb = cbax.colorbar(cs)
cb.set_label('Covariance [C]')

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
cb.set_label('Covariance [C]')

ax.set_title('Model SST / TPI')


########################################################################### ChlSat

'''

# processing obs
data = xr.open_dataset("../chl-sat/interp/covariance_satellite_data.nc")
cov = data['cov'].to_masked_array()
cov = np.ma.masked_where(cov == 0, cov)
lon = data['lon'].values
lat = data['lat'].values

iiii = 4
ax = axgr[iiii]
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

ax.set_title('OCCI Chl / ONI')
cb = axgr.cbar_axes[iiii].colorbar(cs)
cb.set_label("Covariance [mg/m3]")
ax.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

#######

mesh = xr.open_dataset('../data/mesh_mask_eORCA1_v2.2.nc')
mesh = mesh.isel(t=0)
lon = mesh['glamt'].values
lat = mesh['gphit'].values
tmask = mesh['tmask'].values[0]
lonf = mesh['glamf'].values
latf = mesh['gphit'].values

data = xr.open_dataset("../chl-sat/model/covariance_model_data.nc")
cov = data['cov'].values
cov = np.ma.masked_where(tmask == 0, cov)

iiii = 5
ax2 = axgr[iiii]

cs = ax2.pcolormesh(lonf, latf, cov[1:, 1:], transform=proj2)
cs.set_clim(-ccc, ccc)
ax2.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax2.add_feature(cfeature.COASTLINE, zorder=1001)

ax2.set_title('Model Chl / ONI')

xmin = 0.2
#cax = plt.axes([xmin, 0.1, 1-2*xmin, 0.03])
cb = axgr.cbar_axes[iiii].colorbar(cs)
cb.set_label("Covariance [mg/m3]")
ax2.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

gl = ax2.gridlines(**gridparams)
gl.xlabels_top = False
gl.ylabels_right = False
#gl.ylabels_left = False
#gl.xlines = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])

'''


############################################ plotting data


dnino, nino = read_index('../data/index/oni.data')
iok = np.nonzero((dnino >= 195801) & (dnino <= 201812))[0]
dnino = dnino[iok]
nino = nino[iok]
time = np.arange(len(dnino))
ynino = dnino // 100
mnino = dnino - 100 * ynino

labels = np.array(['%.4d-%.2d' %(y, m) for y,m in zip(ynino, mnino)])

data = xr.open_dataset("../nemo-pisces/nino/simulated_enso_index.nc")
years = data['time'].values // 100
months = data['time'].values - 100 * years
enso = data['enso'].values
#timemod = np.arange(len(years)) + 8 * 12
nmod = len(time)
#enso = data['enso'].values
#clim, enso = ts.get_monthly_clim(enso)
#enso = sig.detrend(enso)
ensof = np.zeros(enso.shape)

index = np.arange(3)
for i in range(1, nmod - 1):
    ensof[i] = enso[index].mean()
    index += 1

ensof[ensof == 0] = np.nan

istart = np.nonzero(dnino == 195802)[0][0]
iend = np.nonzero(dnino == 201811)[0][0] + 1

test = np.corrcoef(ensof[1:-1], nino[istart:iend])
print(ensof.shape, nino.shape)
test = np.corrcoef(ensof[1:-1], nino[1:-1])[0, 1]
print('Correlation ONI', test)

istart = np.nonzero(dnino == 195801)[0][0]
iend = np.nonzero(dnino == 201812)[0][0] + 1
test = np.corrcoef(enso, nino[istart:iend])

xticks = np.arange(2*12, len(time), 5 * 12)

left = 0.05
width = 0.4
bottom = 0.6
height = 0.1

axes = (left, bottom, width, height)

alpha = 0.7

ax = plt.axes(axes)
l1 = plt.fill_between(time, 0, nino, where=(nino>0), color='firebrick', interpolate=True)
l2 = plt.fill_between(time, 0, nino, where=(nino<0), color='steelblue', interpolate=True)
l3 = plt.plot(time, ensof, 'k', label='Sim.', alpha=alpha)
plt.legend([l1, l3[0]], ['Obs.', 'Model'], loc=0, fontsize=8, ncol=2)
ax.set_title('ONI')

#plt.legend(loc=0)
ax.set_xticks(time[xticks])
ax.set_xticklabels(labels[xticks], rotation=45, ha='right')
ax.grid(True)
ax.set_xlim(time.min(), time.max())
ax.set_ylim(-4, 4)
ax.text(time[-1] - 50, -3, 'a' + ")", ha='center', va='center', bbox=dicttext)

data = xr.open_dataset('../data/filt_tpi.nc')
nino2 = data['tpi_filt'].values

data = xr.open_dataset('../nemo-pisces/nino/model_tpi_index.nc')
modtpi = data['tpi'].values

modtpi = lanc.wgt_runave(modtpi)
#nino2 = lanc.wgt_runave(nino2)

print('Correlation TPI', np.corrcoef(modtpi, nino2)[0, 1])

left = 0.55
width = 0.4
bottom = 0.6
height = 0.1
axes = (left, bottom, width, height)

ax = plt.axes(axes)
l1 = plt.fill_between(time, 0, nino2, where=(nino2 > 0), color='firebrick', interpolate=True)
l2 = plt.fill_between(time, 0, nino2, where=(nino2 < 0), color='steelblue', interpolate=True)
l3 = plt.plot(time, modtpi, color='k', label='Sim.', alpha=alpha)
plt.legend([l1, l3[0]], ['Obs', 'Model'], loc=0, fontsize=8, ncol=2)
ax.set_title('TPI')

#plt.legend(loc=0)
ax.set_xticks(time[xticks])
ax.set_xticklabels(labels[xticks], rotation=45, ha='right')
ax.grid(True)
ax.set_xlim(time.min(), time.max())
ax.set_ylim(-1, 1)
ax.text(time[-1] - 50, -0.75, 'b' + ")", ha='center', va='center', bbox=dicttext)

####################################################### equatorial tempratures

"""

width = 0.8
left = 0.1
bottom = 0.1
height = 0.12
axes = (left, bottom, width, height)

data = xr.open_dataset('../cov-hadley/obs_equatorial_mean.nc')
obs = data['sst'].values

data = xr.open_dataset('../nemo-pisces/nino/simulated_equatorial_mean.nc')
mod = data['sst'].values

ax = plt.axes(axes)
l1 = plt.plot(time, obs, color='k')
l3 = plt.plot(time, mod, color='firebrick', label='Sim.', alpha=0.5)
plt.legend([l1[0], l3[0]], ['Obs', 'Model'], loc=0, fontsize=8, ncol=2)
ax.set_title('Equatorial mean SST')
ax.set_xticks(time[xticks])
ax.set_xticklabels(labels[xticks], rotation=45, ha='right')
ax.grid(True)
ax.set_xlim(time.min(), time.max())
ax.set_ylabel('[C]')
#ax.set_ylim(-1, 1)
#ax.text(time[-1] - 50, -0.75, 'b' + ")", ha='center', va='center', bbox=dicttext)

"""

plt.savefig('fig1', bbox_inches='tight')
