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
plt.rcParams['font.size'] = 15

latbox = [-5, -5, 5, 5, -5]
lonbox = [-170, -120, -120, -170, -170]


letters = list(string.ascii_lowercase)
letters = letters[1:]

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

dictpbox = {'transform': proj2, 'linestyle': '--', 'linewidth': 1, 'color':'k'}

gridparams = {'crs': ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}

mesh = xr.open_dataset('../../data/mesh_mask_eORCA1_v2.2.nc')
mesh = mesh.isel(t=0, z=0)
lon = mesh['glamt'].values
lat = mesh['gphit'].values
lonf = mesh['glamf'].values
latf = mesh['gphif'].values
tmask = mesh['tmask'].values

ilat = slice(None, -3)

# fig = plt.figure(figsize=(8, 16), facecolor='white')
# plt.subplots_adjust(top=0.95)
# axes_class = (GeoAxes, dict(map_projection=proj))

# left = 0.06
# width = 1 - 2 * left
# bottom = 0.1
# height = 0.45

# axes = (left, bottom, width, height)
# axgr = AxesGrid(fig, axes,  axes_class=axes_class, nrows_ncols=(2, 2), axes_pad=(0.7, 0.75), label_mode='', cbar_mode='each', cbar_size=0.1, cbar_pad=0.3, cbar_location="bottom")

# axout = list(enumerate(axgr))
# axout = [p[1] for p in axout]

# ################################################################### Plotting covariances OBS

# hadley = xr.open_dataset('data/cov_hadsst_oni_tpi.nc')
# lonhad = hadley['lon'].values
# lathad = hadley['lat'].values
# hadoni = hadley['covoni'].to_masked_array().T
# hadoni = np.ma.masked_where(hadoni == 0, hadoni)

# model = xr.open_dataset('data/cov_modssh_oni_tpi.nc')
# modsshoni = model['covoni'].to_masked_array().T * 100
# modsshoni = np.ma.masked_where(modsshoni == 0, modsshoni)

# model = xr.open_dataset('data/cov_modsst_oni_tpi.nc')
# modoni = model['covoni'].to_masked_array().T
# modoni = np.ma.masked_where(modoni == 0, modoni)

# model = xr.open_dataset('data/cov_obsssh_oni_tpi.nc')
# lonsat = model['longitude'].values
# latsat = model['latitude'].values
# obssshoni = model['covoni'].to_masked_array().T * 100
# obssshoni = np.ma.masked_where(obssshoni == 0, obssshoni)

# lontext = 120
# lattext = 30
# dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)

# iiii = 0
# ax = axgr[iiii]
# ccc = 1.1
# step = 0.1
# levels = np.arange(-ccc, ccc + step, step)

# cs = ax.pcolormesh(lonhad, lathad, hadoni, transform=proj2, shading='auto')
# cs.set_clim(-ccc, ccc)
# ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
# ax.add_feature(cfeature.COASTLINE, zorder=1001)
# ax.plot(lonbox, latbox, **dictpbox)

# gl = ax.gridlines(**gridparams)
# gl.top_labels = False
# gl.right_labels = False
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER
# gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])
# ax.set_ylim(-40, 40)
# ax.set_xlim(-60, 130)
# ax.text(lontext, lattext, 'b' + ")", ha='center', va='center', transform=proj, bbox=dicttext)

# ax.set_title('Hadley SST / ONI')
# cbax = axgr.cbar_axes[iiii]
# cb = cbax.colorbar(cs)
# try:
#     cb.set_label('Covariance [C]')
# except:
#     cb.set_label_text('Covariance [C]')
    
# ################################################################### Plotting covariances SST MOD
    
# iiii = 2
# ax = axgr[iiii]

# cs = ax.pcolormesh(lonf, latf, modoni[1:, 1:], transform=proj2)
# cs.set_clim(-ccc, ccc)
# ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
# ax.add_feature(cfeature.COASTLINE, zorder=1001)
# ax.text(lontext, lattext, 'c' + ")", ha='center', va='center', transform=proj, bbox=dicttext)
# ax.plot(lonbox, latbox, **dictpbox)

# gl = ax.gridlines(**gridparams)
# gl.top_labels = False
# gl.right_labels = False
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER
# xticks = np.array([150, 180, -180, -150, -120, -90, -60])
# gl.xlocator = mticker.FixedLocator(xticks)
# ax.set_ylim(-40, 40)
# ax.set_xlim(-60, 130)

# cbax = axgr.cbar_axes[iiii]
# cb = cbax.colorbar(cs)
# try:
#     cb.set_label('Covariance [C]')
# except:
#     cb.set_label_text('Covariance [C]')

# ax.set_title('Model SST / ONI')

# ################################################################### Plotting covariances SSH MOD
# ilat = slice(None, -3)
# modcovv = xr.open_dataset('data/model_covariance_nino34_vo.nc').isel(y=ilat)
# modcovv = modcovv['__xarray_dataarray_variable__']

# modcovu = xr.open_dataset('data/model_covariance_nino34_uo.nc').isel(y=ilat)
# modcovu = modcovu['__xarray_dataarray_variable__']
# modcovu

# modmesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0, y=ilat)
# lonvec = modmesh['glamt'].values
# latvec = modmesh['gphit'].values

# quivargs = {}
# quivargs['scale'] = 1 
# quivargs['width'] = 0.004
# xkey = 0.92
# ykey = 0.32

# iiii = 3
# ax = axgr[iiii]

# cs = ax.pcolormesh(lonf, latf, modsshoni[1:, 1:].T, transform=proj2)
# q = ax.quiver(lonvec, latvec, modcovu.values, modcovv.values, transform=ccrs.PlateCarree(), **quivargs, regrid_shape=(30, 30))
# ref = 0.05
# plt.quiverkey(q, xkey, ykey , ref, '%dcm/s' %(ref * 100) , color='k')

# #cs.set_clim(-ccc, ccc)
# ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
# ax.add_feature(cfeature.COASTLINE, zorder=1001)
# ax.text(lontext, lattext, 'f' + ")", ha='center', va='center', transform=proj, bbox=dicttext)
# ax.plot(lonbox, latbox, **dictpbox)

# gl = ax.gridlines(**gridparams)
# gl.top_labels = False
# gl.right_labels = False
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER
# xticks = np.array([150, 180, -180, -150, -120, -90, -60])
# gl.xlocator = mticker.FixedLocator(xticks)
# ax.set_ylim(-40, 40)
# ax.set_xlim(-60, 130)
# cs.set_clim(-8, 8)


# cbax = axgr.cbar_axes[iiii]
# cb = cbax.colorbar(cs)
# try:
#     cb.set_label('Covariance [cm]')
# except:
#     cb.set_label_text('Covariance [cm]')

# ax.set_title('Model SSH / ONI')

# ################################################################### Plotting covariances SSH OBS

# satcovu = xr.open_dataset('data/sat_covariance_nino34_uo.nc')
# satcovu = satcovu['__xarray_dataarray_variable__']
# satcovu

# satcovv = xr.open_dataset('data/sat_covariance_nino34_vo.nc')
# satcovv = satcovv['__xarray_dataarray_variable__']
# satcovv

# lonobs = satcovu['longitude'].values
# latobs = satcovu['latitude'].values
    
# iiii = 1
# ax = axgr[iiii]

# cs = ax.pcolormesh(lonsat, latsat, obssshoni[:, :].T, transform=proj2, shading='auto')
# q = ax.quiver(lonobs, latobs, satcovu.values, satcovv.values, transform=ccrs.PlateCarree(), **quivargs, regrid_shape=(20, 20))
# ref = 0.05
# plt.quiverkey(q, xkey, ykey , ref, '%dcm/s' %(ref * 100) , color='k')

# ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
# ax.add_feature(cfeature.COASTLINE, zorder=1001)
# ax.text(lontext, lattext, 'e' + ")", ha='center', va='center', transform=proj, bbox=dicttext)
# ax.plot(lonbox, latbox, **dictpbox)

# gl = ax.gridlines(**gridparams)
# gl.top_labels = False
# gl.right_labels = False
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER
# xticks = np.array([150, 180, -180, -150, -120, -90, -60])
# gl.xlocator = mticker.FixedLocator(xticks)
# ax.set_ylim(-40, 40)
# ax.set_xlim(-60, 130)
# cs.set_clim(-8, 8)

# cbax = axgr.cbar_axes[iiii]
# cb = cbax.colorbar(cs)
# try:
#     cb.set_label('Covariance [cm]')
# except:
#     cb.set_label_text('Covariance [cm]')

# ax.set_title('Obs SSH / ONI')


# ################################################################### Plotting time-series ONI index

# dnino, nino = read_index('../../data/index/oni.data')
# iok = np.nonzero((dnino >= 195801) & (dnino <= 201812))[0]
# dnino = dnino[iok]
# nino = nino[iok]
# time = np.arange(len(dnino))
# ynino = dnino // 100
# mnino = dnino - 100 * ynino

# labels = np.array(['%.4d-%.2d' %(y, m) for y,m in zip(ynino, mnino)])

# data = xr.open_dataset("data/simulated_enso_index.nc")
# years = data['time_counter.year'].values 
# months = data['time_counter.month'].values
# enso = data['enso'].rolling(time_counter=3, center=True).mean()
# nmod = len(time)

# istart = np.nonzero(dnino == 195802)[0][0]
# iend = np.nonzero(dnino == 201811)[0][0] + 1

# test = np.corrcoef(enso[1:-1], nino[istart:iend])[0, 1]
# print('Correlation ONI', test)

# istart = np.nonzero(dnino == 195801)[0][0]
# iend = np.nonzero(dnino == 201812)[0][0] + 1
# test = np.corrcoef(enso, nino[istart:iend])

# xticks = np.arange(2*12, len(time), 5 * 12)

# left = 0.05
# width = 0.4
# bottom = 0.53
# height = 0.1

# axes = (left, bottom, width, height)

# alpha = 0.7

# ax = plt.axes(axes)
# l1 = plt.fill_between(time, 0, nino, where=(nino>0), color='firebrick', interpolate=True)
# l2 = plt.fill_between(time, 0, nino, where=(nino<0), color='steelblue', interpolate=True)
# l3 = plt.plot(time, enso, 'k', label='Sim.', alpha=alpha)
# plt.legend([l1, l3[0]], ['Obs.', 'Model'], loc=0, fontsize=8, ncol=2)
# ax.set_title('ONI')

# #plt.legend(loc=0)
# ax.set_xticks(time[xticks])
# ax.set_xticklabels(labels[xticks], rotation=45, ha='right')
# ax.grid(True)
# ax.set_xlim(time.min(), time.max())
# ax.set_ylim(-4, 4)
# ax.text(time[-1] - 50, -3, 'a' + ")", ha='center', va='center', bbox=dicttext)

# ####################################### Plotting SSH time-series

# curobs = xr.open_dataset('data/satellite_nino_34_uo.nc')
# uobs = curobs['uo']
# uobs_clim = uobs.groupby('time.month').mean('time')
# uobs = uobs.groupby('time.month') - uobs_clim

# curobs = xr.open_dataset('data/satellite_nino_34_vo.nc')
# vobs = curobs['vo']
# vobs_clim = uobs.groupby('time.month').mean('time')
# vobs = uobs.groupby('time.month') - vobs_clim

# curmod = xr.open_dataset('data/model_nino_34_uo.nc')
# umod = curmod['__xarray_dataarray_variable__']  
# umod_clim = umod.groupby('time_counter.month').mean('time_counter')
# umod = umod.groupby('time_counter.month') - umod_clim

# curmod = xr.open_dataset('data/model_nino_34_vo.nc')
# vmod = curmod['__xarray_dataarray_variable__']
# vmod_clim = umod.groupby('time_counter.month').mean('time_counter')
# vmod = umod.groupby('time_counter.month') - vmod_clim

# dateobs = np.array([(y * 100 + m) for y, m in zip(uobs['time.year'].values, uobs['time.month'].values)])
# datemod = np.array([(y * 100 + m) for y, m in zip(umod['time_counter.year'].values, umod['time_counter.month'].values)])
# tmod = np.arange(len(umod))

# offset = np.nonzero(datemod == 199301)[0][0]
# tobs = np.arange(len(uobs)) + offset

# corrumod = umod.sel(time_counter=slice('1993-01-01', '2018-12-31'))
# corrumod

# corruobs = uobs.sel(time=slice('1993-01-01', '2018-12-31'))
# corruobs

# print('Correlation U', np.corrcoef(corruobs.values, corrumod.values)[0, 1])

# left = 0.55

# axes = (left, bottom, width, height)

# alpha = 0.7

# ax = plt.axes(axes)
# l3 = plt.plot(tmod, umod, 'k', label='Sim.', alpha=alpha)
# l1 = plt.plot(tobs, uobs, label='Sim.', alpha=alpha, color='orange')
# plt.legend([l1[0], l3[0]], ['Obs.', 'Model'], loc=0, fontsize=8, ncol=2)
# plt.ylabel('[m/s]')
# plt.title('U Nino34')

# #plt.legend(loc=0)
# ax.set_xticks(time[xticks])
# ax.set_xticklabels(labels[xticks], rotation=45, ha='right')
# ax.grid(True)
# ax.set_xlim(time.min(), time.max())
# ccc = 0.7
# ax.set_ylim(-ccc, ccc)
# ax.text(time[-1] - 50, -0.5, 'd' + ")", ha='center', va='center', bbox=dicttext)


# figplt.savefig('fig1', bbox_inches='tight')
# -
# # SST

# +
fig = plt.figure(facecolor='white', figsize=(12, 8))
ax = plt.gca(projection=proj)

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

model = xr.open_dataset('data/cov_obsssh_oni_tpi.nc')
lonsat = model['longitude'].values
latsat = model['latitude'].values
obssshoni = model['covoni'].to_masked_array().T * 100
obssshoni = np.ma.masked_where(obssshoni == 0, obssshoni)

lontext = 120
lattext = 30
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)

iiii = 0
ccc = 1.1
step = 0.1
levels = np.arange(-ccc, ccc + step, step)

cs = ax.pcolormesh(lonhad, lathad, hadoni, transform=proj2, shading='auto')
cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax.add_feature(cfeature.COASTLINE, zorder=1001)
ax.plot(lonbox, latbox, **dictpbox)

gl = ax.gridlines(**gridparams)
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90, -60])
ax.set_ylim(-40, 40)
ax.set_xlim(-60, 130)

ax.set_title('Hadley SST / ONI')
cb = plt.colorbar(cs, shrink=0.5)
try:
    cb.set_label('Covariance [C]')
except:
    cb.set_label_text('Covariance [C]')
plt.savefig('obs_sst_cov.png', bbox_inches='tight')

# +
fig = plt.figure(facecolor='white', figsize=(12, 8))
ax = plt.gca(projection=proj)

cs = ax.pcolormesh(lonf, latf, modoni[1:, 1:], transform=proj2)
cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
ax.add_feature(cfeature.COASTLINE, zorder=1001)
ax.plot(lonbox, latbox, **dictpbox)

gl = ax.gridlines(**gridparams)
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
xticks = np.array([150, 180, -180, -150, -120, -90, -60])
gl.xlocator = mticker.FixedLocator(xticks)
ax.set_ylim(-40, 40)
ax.set_xlim(-60, 130)

cb = plt.colorbar(cs, shrink=0.5)
try:
    cb.set_label('Covariance [C]')
except:
    cb.set_label_text('Covariance [C]')

ax.set_title('Model SST / ONI')
plt.savefig('mod_sst_cov.png', bbox_inches='tight')

# +
plt.figure(facecolor='white', figsize=(14, 8))
ax = plt.gca()

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
nmod = len(time)

istart = np.nonzero(dnino == 195802)[0][0]
iend = np.nonzero(dnino == 201811)[0][0] + 1

test = np.corrcoef(enso[1:-1], nino[istart:iend])[0, 1]
print('Correlation ONI', test)

istart = np.nonzero(dnino == 195801)[0][0]
iend = np.nonzero(dnino == 201812)[0][0] + 1
test = np.corrcoef(enso, nino[istart:iend])

xticks = np.arange(2*12, len(time), 5 * 12)

alpha = 0.7

l1 = plt.fill_between(time, 0, nino, where=(nino>0), color='firebrick', interpolate=True)
l2 = plt.fill_between(time, 0, nino, where=(nino<0), color='steelblue', interpolate=True)
l3 = plt.plot(time, enso, 'k', label='Sim.', alpha=1)
plt.legend([l1, l3[0]], ['Obs.', 'Model'], loc=0, ncol=2)
ax.set_title('ENSO index')

#plt.legend(loc=0)
ax.set_xticks(time[xticks])
ax.set_xticklabels(labels[xticks], rotation=45, ha='right')
ax.grid(True)
ax.set_xlim(time.min(), time.max())
ax.set_ylim(-4, 4)
#ax.text(time[-1] - 50, -3, 'a' + ")", ha='center', va='center', bbox=dicttext)
plt.savefig('oni_index.png', bbox_inches='tight')
# -
# # SSH / Currents

# +
quivargs = {}
quivargs['scale'] = 1
quivargs['width'] = 0.004
xkey = 0.92
ykey = 0.32

satcovu = xr.open_dataset('data/sat_covariance_nino34_uo.nc')
satcovu = satcovu['__xarray_dataarray_variable__']
satcovu

satcovv = xr.open_dataset('data/sat_covariance_nino34_vo.nc')
satcovv = satcovv['__xarray_dataarray_variable__']
satcovv

lonobs = satcovu['longitude'].values
latobs = satcovu['latitude'].values
    
iiii = 1
fig = plt.figure(facecolor='white', figsize=(12, 8))
ax = plt.axes(projection=proj)

cs = ax.pcolormesh(lonsat, latsat, obssshoni[:, :].T, transform=proj2, shading='auto')
q = ax.quiver(lonobs, latobs, satcovu.values, satcovv.values, transform=ccrs.PlateCarree(), **quivargs, regrid_shape=(20, 20))
ref = 0.05

ax.add_feature(cfeature.LAND, color='lightgray')
ax.add_feature(cfeature.COASTLINE)
ax.text(lontext, lattext, 'e' + ")", ha='center', va='center', transform=proj, bbox=dicttext)
ax.plot(lonbox, latbox, **dictpbox)

ax.quiverkey(q, xkey, ykey , ref, '%dcm/s' %(ref * 100) , color='k', zorder=1003)

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

cb = plt.colorbar(cs, shrink=0.5)
try:
    cb.set_label('Covariance [cm]')
except:
    cb.set_label_text('Covariance [cm]')

ax.set_title('Obs SSH / ONI')

plt.savefig('obs_ssh_uv_cov.png', bbox_inches='tight')

# +
modcovv = xr.open_dataset('data/model_covariance_nino34_vo.nc').isel(y=ilat)
modcovv = modcovv['__xarray_dataarray_variable__']

modcovu = xr.open_dataset('data/model_covariance_nino34_uo.nc').isel(y=ilat)
modcovu = modcovu['__xarray_dataarray_variable__']
modcovu

modmesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0, y=ilat)
lonvec = modmesh['glamt'].values
latvec = modmesh['gphit'].values

iiii = 3
fig = plt.figure(facecolor='white', figsize=(12, 8))
ax = plt.axes(projection=proj)

cs = ax.pcolormesh(lonf, latf, modsshoni[1:, 1:].T, transform=proj2)
q = ax.quiver(lonvec, latvec, modcovu.values, modcovv.values, transform=ccrs.PlateCarree(), **quivargs, regrid_shape=(40, 40))
ref = 0.05

#cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, color='lightgray')
ax.add_feature(cfeature.COASTLINE)
ax.text(lontext, lattext, 'f' + ")", ha='center', va='center', transform=proj, bbox=dicttext)
ax.plot(lonbox, latbox, **dictpbox)

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

plt.quiverkey(q, xkey, ykey , ref, '%dcm/s' %(ref * 100) , color='k', zorder=30000)

cb = plt.colorbar(cs, shrink=0.5)
try:
    cb.set_label('Covariance [cm]')
except:
    cb.set_label_text('Covariance [cm]')

ax.set_title('Model SSH / ONI')
plt.savefig('mod_ssh_uv_cov.png', bbox_inches='tight')


# +
curobs = xr.open_dataset('data/satellite_nino_34_uo.nc')
uobs = curobs['uo']
uobs_clim = uobs.groupby('time.month').mean('time')
uobs = uobs.groupby('time.month') - uobs_clim

curobs = xr.open_dataset('data/satellite_nino_34_vo.nc')
vobs = curobs['vo']
vobs_clim = uobs.groupby('time.month').mean('time')
vobs = uobs.groupby('time.month') - vobs_clim

curmod = xr.open_dataset('data/model_nino_34_uo.nc')
umod = curmod['__xarray_dataarray_variable__']  
umod_clim = umod.groupby('time_counter.month').mean('time_counter')
umod = umod.groupby('time_counter.month') - umod_clim

curmod = xr.open_dataset('data/model_nino_34_vo.nc')
vmod = curmod['__xarray_dataarray_variable__']
vmod_clim = umod.groupby('time_counter.month').mean('time_counter')
vmod = umod.groupby('time_counter.month') - vmod_clim

dateobs = np.array([(y * 100 + m) for y, m in zip(uobs['time.year'].values, uobs['time.month'].values)])
datemod = np.array([(y * 100 + m) for y, m in zip(umod['time_counter.year'].values, umod['time_counter.month'].values)])
tmod = np.arange(len(umod))

offset = np.nonzero(datemod == 199301)[0][0]
tobs = np.arange(len(uobs)) + offset

corrumod = umod.sel(time_counter=slice('1993-01-01', '2018-12-31'))
corrumod

corruobs = uobs.sel(time=slice('1993-01-01', '2018-12-31'))
corruobs

print('Correlation U', np.corrcoef(corruobs.values, corrumod.values)[0, 1])

left = 0.55

alpha = 0.7

xticks = np.arange(2*12, len(time), 5 * 12)

plt.figure(facecolor='white', figsize=(14, 8))
ax = plt.gca()
l3 = plt.plot(tmod, umod, 'k', label='Sim.', alpha=alpha)
l1 = plt.plot(tobs, uobs, label='Sim.', alpha=alpha, color='orange')
plt.legend([l1[0], l3[0]], ['Obs.', 'Model'], loc=0, fontsize=15, ncol=2)
plt.ylabel('[m/s]')
plt.title('U Nino34')

#plt.legend(loc=0)
ax.set_xticks(time[xticks])
ax.set_xticklabels(labels[xticks], rotation=45, ha='right')
ax.grid(True)
ax.set_xlim(time.min(), time.max())
ccc = 0.7
ax.set_ylim(-0.4, ccc)
ax.text(time[-1] - 50, -0.5, 'd' + ")", ha='center', va='center', bbox=dicttext)

plt.savefig('u_nino34.png', bbox_inches='tight')
# -


