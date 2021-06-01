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
from envtoolkit.ts import Lanczos

nWgt      = 157
fca       = 1./156
fca = 1 / fca

lanc = Lanczos('lp', pca=fca, nwts=nWgt)

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)

letters = list(string.ascii_lowercase)
letters = letters[2:]

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

fig = plt.figure(figsize=(8, 16))

plt.subplots_adjust(top=0.95)
axes_class = (GeoAxes, dict(map_projection=proj))

left = 0.06
width = 1 - 2 * left
bottom = 0.1
height = 0.45

axes = (left, bottom, width, height)
axgr = AxesGrid(fig, axes,  axes_class=axes_class, nrows_ncols=(1, 2), axes_pad=(0.7, 0.75), label_mode='', cbar_mode='each', cbar_size=0.1, cbar_pad=0.3, cbar_location="bottom")

axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

#################################################################### Plotting covariances satellite

ccc = 0.1
lontext = 120
lattext = 30

# processing obs
data = xr.open_dataset("../../chl-sat/interp/covariance_satellite_data.nc")
cov = data['cov'].to_masked_array()
cov = np.ma.masked_where(cov == 0, cov)
lon = data['lon'].values
lat = data['lat'].values

iiii = 0
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
cb.set_label_text("Covariance [mg/m3]")
ax.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

###################################### plot cov. model

mesh = xr.open_dataset('../../data/mesh_mask_eORCA1_v2.2.nc')
mesh = mesh.isel(t=0)
lon = mesh['glamt'].values
lat = mesh['gphit'].values
tmask = mesh['tmask'].values[0]
lonf = mesh['glamf'].values
latf = mesh['gphit'].values

data = xr.open_dataset("../../chl-sat/model/covariance_model_data.nc")
cov = data['cov'].values
cov = np.ma.masked_where(tmask == 0, cov)

data = xr.open_dataset("../../chl-sat/model/cov_modchl_oni_tpi.nc")
cov = data['covoni'].values
cov = cov.T
cov = np.ma.masked_where(tmask == 0, cov)

iiii = 1
ax2 = axgr[iiii]

cs = ax2.pcolormesh(lonf, latf, cov[1:, 1:], transform=proj2)
cs.set_clim(-ccc, ccc)
ax2.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax2.add_feature(cfeature.COASTLINE, zorder=1001)
ax2.set_ylim(-40, 40)
ax2.set_xlim(-60, 130)

ax2.set_title('Model Chl / ONI')

xmin = 0.2
#cax = plt.axes([xmin, 0.1, 1-2*xmin, 0.03])
cb = axgr.cbar_axes[iiii].colorbar(cs)
cb.set_label_text("Covariance [mg/m3]")
ax2.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

gl = ax2.gridlines(**gridparams)
gl.xlabels_top = False
gl.ylabels_right = False
#gl.ylabels_left = False
#gl.xlines = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])

################################################ tpi
left = 0.54
width = 0.4
bottom = 0.09
height = 0.15
ccc2 = 1e-2

axes = (left, bottom, width, height)
data = xr.open_dataset("../../chl-sat/model/cov_modchl_oni_tpi.nc")
cov = data['covtpi'].values
cov = cov.T
cov = np.ma.masked_where(tmask == 0, cov)

iiii = 1
ax2 = plt.axes(axes, projection=proj)

cs = ax2.pcolormesh(lonf, latf, cov[1:, 1:], transform=proj2)
cs.set_clim(-ccc2, ccc2)
ax2.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax2.add_feature(cfeature.COASTLINE, zorder=1001)
ax2.set_ylim(-40, 40)
ax2.set_xlim(-60, 130)

ax2.set_title('Model Chl / TPI')

xmin = 0.2
#cax = plt.axes([xmin, 0.1, 1-2*xmin, 0.03])
cb = plt.colorbar(cs, orientation='horizontal', shrink=0.8)
cb.set_label("Covariance [mg/m3]")
ax2.text(lontext, lattext, letters[2] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

gl = ax2.gridlines(**gridparams)
gl.xlabels_top = False
gl.ylabels_right = False
#gl.ylabels_left = False
#gl.xlines = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])

############################################ plotting time-series

data = xr.open_dataset('../../data/filt_tpi.nc')
tpi = data['tpi_filt'].values

dnino, nino = read_index('../../data/index/oni.data')
iok = np.nonzero((dnino >= 195801) & (dnino <= 201812))[0]
dnino = dnino[iok]
nino = nino[iok]

data = xr.open_dataset('../../chl-sat/model/simulated_equatorial_mean.nc')
datemod = data['time_counter.year'] * 100 + data['time_counter.month'] 
datemod = datemod.values
mod = data['chl'].values

data = xr.open_dataset('../../chl-sat/interp/obs_equatorial_mean.nc')
obs = data['chl'].values
date = data['time.year'] * 100 + data['time.month']
date = date.values
otime = np.arange(len(mod))

imod = np.nonzero(datemod >= date[0])[0]
iobs = np.nonzero(date <= datemod[-1])[0]

toffset = imod[0]

time = np.arange(len(obs))
timemod = np.arange(len(mod))

print(date)
iclimobs = np.nonzero((date >= 199801) & (date <= 201812))[0]
iclimmod = np.nonzero((datemod >= 195801) & (datemod <= 201812))[0]
climmod, anom = ts.get_monthly_clim(mod[iclimmod])
climobs, anom = ts.get_monthly_clim(obs[iclimobs])
year  = date // 100
month = date - 100 * year

yearmod  = datemod // 100
monthmod = datemod - 100 * yearmod

labels = ['%.4d-%.2d' %(y, m) for y, m in zip(yearmod, monthmod)]
labels = np.array(labels)

for t in time:
    m = month[t] - 1
    obs[t] = obs[t] - climobs[m] 

for t in timemod:
    m = monthmod[t] - 1
    mod[t] = mod[t] - climmod[m] 

print(np.mean(mod))

print('Correlation = ', np.corrcoef(obs[iobs], mod[imod])[0, 1])

xticks = np.arange(2*12, len(timemod), 5*12)

left = 0.07
width = 0.37
bottom = 0.45
height = 0.1

axes = (left, bottom, width, height)

alpha = 0.7

ax = plt.axes(axes)
l2 = plt.plot(timemod, mod, 'black', alpha=alpha)
l1 = plt.plot(time + toffset, obs, color='orange', alpha=alpha)
ax.set_xlim(toffset, timemod.max())
leg = plt.legend([l1[0], l2[0]], ['Obs.', 'Model'], loc=0, fontsize=8, ncol=2)
ax.add_artist(leg)
ax.set_title('Equatorial CHLA anomalies')
ax.set_ylabel('[mg/m3]')

#plt.legend(loc=0)
ax.set_xticks(timemod[xticks])
ax.set_xticklabels(labels[xticks], rotation=45, ha='right')
ax.grid(True)
#ax.set_xlim(time.min(), time.max())
ax.set_ylim(-0.2, 0.2)
ax.text(timemod[-1] - 50, -0.15, 'a' + ")", ha='center', va='center', bbox=dicttext)

#axbis = ax.twinx()
#ax.set_zorder(1)
#ax.patch.set_visible(False)
#lll = axbis.fill_between(timemod, 0, nino, where=(nino > 0), color='firebrick', interpolate=True)
#lll = axbis.fill_between(timemod, 0, nino, where=(nino < 0), color='steelblue', interpolate=True)
#mmm = 2.5
#axbis.set_ylim(-mmm, mmm)
#axbis.get_yaxis().set_visible(False)  # removes xlabels

left = 0.56

axes = (left, bottom, width, height)

modf = lanc.wgt_runave(mod)

alpha = 0.7

ax = plt.axes(axes)
l2 = plt.plot(timemod, modf, 'black', alpha=alpha)
ax.set_xlim(toffset, timemod.max())
leg = plt.legend([l2[0]], ['Model'], loc=0, fontsize=8, ncol=2)
ax.add_artist(leg)
ax.set_title('Equatorial CHLA anomalies')
ax.set_ylabel('[mg/m3]')

ax.set_xticks(timemod[xticks])
ax.set_xticklabels(labels[xticks], rotation=45, ha='right')
ax.grid(True)
ax.set_xlim(timemod.min(), timemod.max())
yyy = 0.028
ax.set_ylim(-yyy, yyy)
ax.text(timemod[-1] - 50, -0.0175, 'b' + ")", ha='center', va='center', bbox=dicttext)

#axbis = ax.twinx()
#ax.set_zorder(2)
#ax.patch.set_visible(False)
#lll = axbis.fill_between(timemod, 0, tpi, where=(tpi > 0), color='firebrick', interpolate=True)
#lll = axbis.fill_between(timemod, 0, tpi, where=(tpi < 0), color='steelblue', interpolate=True)
#axbis.set_ylim(-1, 1)
##axbis.set_ylabel('TPI')
##legend = plt.legend([lll[0]], ['TPI'], loc=2, fontsize=8)
#axbis.get_yaxis().set_visible(False)  # removes xlabels

plt.savefig('fig2', bbox_inches='tight')

