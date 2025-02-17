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
sys.path.append('..')
from extract_nino import read_index
import scipy.signal as sig
plt.rcParams['image.cmap'] = 'RdBu_r'

y = slice(None, -3)

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)

letters = list(string.ascii_lowercase)
letters = letters[1:]

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

latbox = [-2, -2, 2, 2, -2]
lonbox = [150, 150, -80, -80, 150]
lonbox = [150, -80 + 360, -80 + 360, 150, 150]
dictpbox = {'transform': proj2, 'linestyle': '--', 'linewidth': 1, 'color':'k'}

gridparams = {'crs': ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}

mesh = xr.open_dataset('../data/static/pacific_mesh_mask.nc').isel(y=y, z=0)
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
axgr = AxesGrid(fig, axes,  axes_class=axes_class, nrows_ncols=(2, 1), axes_pad=(0.7, 0.75), label_mode='', cbar_mode='each', cbar_size=0.1, cbar_pad=0.3, cbar_location="bottom")

axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

################################################################### Plotting covariances satellite

ccc = 0.1
lontext = 120
lattext = 30

# processing obs
data = xr.open_dataset("chl-sat/covariance_satellite_data.nc")
cov = data['cov'].to_masked_array()
cov = np.ma.masked_where(cov == 0, cov)
lon = data['lon'].values
lat = data['lat'].values

iiii = 0
ax = axgr[iiii]
ax.set_ylim(-40, 40)
ax.set_xlim(-60, 130)

cs = ax.pcolormesh(lon, lat, cov, transform=ccrs.PlateCarree(), shading='auto')
ax.plot(lonbox, latbox, **dictpbox)

cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax.add_feature(cfeature.COASTLINE, zorder=1001)

gl = ax.gridlines(**gridparams)
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])

ax.set_title('OCCI Chl / ONI')
cb = axgr.cbar_axes[iiii].colorbar(cs)
cb.set_label("Covariance [mg/m3]")
ax.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

##################################### plot cov. model

# mesh = xr.open_dataset('../../data/mesh_mask_eORCA1_v2.2.nc')
# mesh = mesh.isel(t=0)
# lon = mesh['glamt'].values
# lat = mesh['gphit'].values
# tmask = mesh['tmask'].values[0]
# lonf = mesh['glamf'].values
# latf = mesh['gphit'].values

data = xr.open_dataset("chl-mod/cov_modchl_oni_tpi.nc").isel(y=y)
cov = data['covoni'].values
cov = np.ma.masked_where(tmask == 0, cov)

iiii = 1
ax2 = axgr[iiii]

cs = ax2.pcolormesh(lonf, latf, cov[1:, 1:], transform=proj2)
ax2.plot(lonbox, latbox, **dictpbox)
cs.set_clim(-ccc, ccc)
ax2.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax2.add_feature(cfeature.COASTLINE, zorder=1001)
ax2.set_ylim(-40, 40)
ax2.set_xlim(-60, 130)

ax2.set_title('Model Chl / ONI')

xmin = 0.2
#cax = plt.axes([xmin, 0.1, 1-2*xmin, 0.03])
cb = axgr.cbar_axes[iiii].colorbar(cs)
cb.set_label("Covariance [mg/m3]")
ax2.text(lontext, lattext, letters[iiii] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

gl = ax2.gridlines(**gridparams)
gl.top_labels = False
gl.right_labels = False
#gl.ylabels_left = False
#gl.xlines = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])

########################################### plotting time-series

iclim= slice('1998-01-01', '2018-12-31')

data = xr.open_dataset('chl-mod/simulated_equatorial_mean.nc')
datemod = data['time_counter.year'] * 100 + data['time_counter.month'] 
datemod = datemod.values
mod = data['chl']
mod = mod.groupby('time_counter.month') - mod.sel(time_counter=iclim).groupby('time_counter.month').mean(dim='time_counter')

data = xr.open_dataset('chl-sat/obs_equatorial_mean.nc')
obs = data['chl']
obs = obs.groupby('time.month') - obs.sel(time=iclim).groupby('time.month').mean(dim='time')
date = obs['time.year'] * 100 + obs['time.month']
date = date.values
otime = np.arange(len(mod))

# imod = index of simulation time-steps common with data
imod = np.nonzero(datemod >= date[0])[0]

# iobs = index of data time-steps common with simulation
iobs = np.nonzero(date <= datemod[-1])[0]

timemod = np.arange(len(mod))

yearmod  = datemod // 100
monthmod = datemod - 100 * yearmod

labels = ['%.4d-%.2d' %(y, m) for y, m in zip(yearmod, monthmod)]
labels = np.array(labels)

print('Correlation = ', np.corrcoef(obs[iobs], mod[imod])[0, 1])
ts1 = obs[iobs]
ts2 = mod[imod]
# test = signi.sig(ts1, ts2, dof=14)
# print('Test = ', test)
    
xticks = np.arange(2*12, len(timemod), 5*12)

left = 0.07
width = 0.85
bottom = 0.62
height = 0.1

axes = (left, bottom, width, height)

alpha = 0.7

ax = plt.axes(axes)
l2 = plt.plot(timemod, mod, 'black', alpha=alpha)
l1 = plt.plot(timemod[imod], obs[iobs], color='orange', alpha=alpha)
#ax.set_xlim(toffset, timemod.max())
leg = plt.legend([l1[0], l2[0]], ['Obs.', 'Model'], loc=0, fontsize=8, ncol=2)
ax.add_artist(leg)
ax.set_title('Equatorial CHLA anomalies')
ax.set_ylabel('[mg/m3]')

ax.set_xticks(timemod[xticks])
ax.set_xticklabels(labels[xticks], rotation=45, ha='right')
ax.grid(True)
#ax.set_xlim(time.min(), time.max())
ax.set_ylim(-0.2, 0.2)
ax.text(timemod[-1] -20, -0.15, 'a' + ")", ha='center', va='center', bbox=dicttext)

axbis = ax.twinx()
ax.set_zorder(1)
ax.patch.set_visible(False)
mmm = 2.5
axbis.set_ylim(-mmm, mmm)
axbis.get_yaxis().set_visible(False)  # removes xlabels

left = 0.56

axes = (left, bottom, width, height)

alpha = 0.7

plt.savefig('gr2.jpg', bbox_inches='tight')
# -


