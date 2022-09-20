# ---
# jupyter:
#   jupytext:
#     formats: py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes

l = 90
itime = slice('1997-10-01', '1997-12-31')
chunk = {'x' : 50, 'y' : 50}
# -

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0)
lonf = mesh['glamf'].values
latf = mesh['gphif'].values

const = xr.open_dataset('data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l', 'gpred': 'c'})
length = const['length'] * 100
const['l'] = length
weight_step = const['weight_step']
weight_step

data = xr.open_dataset('data/pacific_nino97_OOPE.nc').chunk(chunk)
oope = data['OOPE']
oope

dataclim = xr.open_dataset('data/pacific_clim_OOPE.nc').chunk(chunk)
oopeclim = dataclim['OOPE']
oopeclim

oopeanom = oope.groupby('time.month') - oopeclim
oopeanom = oopeanom.rename({'w' : 'l'})
oopeanom['l'] = length
oopeanom

oopemean = oopeanom.sel(time=itime).sel(l=l, method='nearest').mean(dim='time')
oopemean

varnames = [
    'growthTrend',
    'madv_trend',
    'mdiff_trend',
    'predationTrend',
    'zadv_trend',
    'zdiff_trend'
]
varnames


def compute_trend(varnames):
    print(varnames)
    for v in varnames:
        temp = xr.open_dataset('data/integrated_maps_%s.nc' %v)
        if(v == varnames[0]):
            trend = temp[v]
        else:
            trend += temp[v]
    trend = trend + oopeanom.isel(time=0)
    trendmean = trend.sel(time=itime).sel(l=l, method='nearest').mean(dim='time')
    trendmean = trendmean.where(trendmean != 0)
    return trendmean


def find_limits(data, thres=95):
    import numpy as np
    temp = np.abs(data[data.mask == False])
    cmax = np.percentile(temp, thres)
    return cmax


# +
projout = ccrs.PlateCarree(central_longitude=180)
projin = ccrs.PlateCarree(central_longitude=0)
axes_class = (GeoAxes, dict(map_projection=projout))

plt.rcParams['font.size'] = 15
plt.rcParams['image.cmap'] = 'RdBu_r'

fig = plt.figure(figsize=(12, 8), facecolor='white')
axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(3, 2), axes_pad=(0.7, 0.95), label_mode='', 
                cbar_mode='each', cbar_size=0.1, cbar_pad=0.3, cbar_location="bottom")

cpt = 0
toplot = oopemean.to_masked_array()[:-3, :][1:, 1:]
cmax = find_limits(toplot)
cs = axgr[cpt].pcolormesh(lonf[:-3, :], latf[:-3, :], toplot, transform=projin)
cbax = axgr.cbar_axes[cpt]
cb = cbax.colorbar(cs)
axgr[cpt].add_feature(cfeature.LAND)
axgr[cpt].add_feature(cfeature.COASTLINE)
cs.set_clim(-cmax, cmax)
axgr[cpt].set_title('OOPE')
cb.set_label('J/m2')

cpt += 1
varnames = [
    'growthTrend',
    'madv_trend',
    'mdiff_trend',
    'predationTrend',
    'zadv_trend',
    'zdiff_trend'
]
toplot = compute_trend(varnames).to_masked_array()[:-3, :][1:, 1:]
cs = axgr[cpt].pcolormesh(lonf[:-3, :], latf[:-3, :], toplot, transform=projin)
cbax = axgr.cbar_axes[cpt]
cb = cbax.colorbar(cs)
axgr[cpt].add_feature(cfeature.LAND)
axgr[cpt].add_feature(cfeature.COASTLINE)
axgr[cpt].set_extent([130, -60, -40, 40], crs=projout)
cmax = find_limits(toplot)
cs.set_clim(-cmax, cmax)
axgr[cpt].set_title('All')
cb.set_label('J/m2')

cpt += 1
varnames = [
    'madv_trend',
    'mdiff_trend',
    'zadv_trend',
    'zdiff_trend'
]
toplot = compute_trend(varnames).to_masked_array()[:-3, :][1:, 1:]
cs = axgr[cpt].pcolormesh(lonf[:-3, :], latf[:-3, :], toplot, transform=projin)
cbax = axgr.cbar_axes[cpt]
cb = cbax.colorbar(cs)
axgr[cpt].add_feature(cfeature.LAND)
axgr[cpt].add_feature(cfeature.COASTLINE)
axgr[cpt].set_extent([130, -60, -40, 40], crs=projout)
cmax = find_limits(toplot)
cs.set_clim(-cmax, cmax)
axgr[cpt].set_title('Movements')
cb.set_label('J/m2')

cpt += 1
varnames = [
    'growthTrend'
]
toplot = compute_trend(varnames).to_masked_array()[:-3, :][1:, 1:]
cs = axgr[cpt].pcolormesh(lonf[:-3, :], latf[:-3, :], toplot, transform=projin)
cbax = axgr.cbar_axes[cpt]
cb = cbax.colorbar(cs)
axgr[cpt].add_feature(cfeature.LAND)
axgr[cpt].add_feature(cfeature.COASTLINE)
axgr[cpt].set_extent([130, -60, -40, 40], crs=projout)
cmax = find_limits(toplot)
cs.set_clim(-cmax, cmax)
axgr[cpt].set_title('Growth')
cb.set_label('J/m2')

cpt += 1
varnames = [
    'predationTrend'
]
toplot = compute_trend(varnames).to_masked_array()[:-3, :][1:, 1:]
cs = axgr[cpt].pcolormesh(lonf[:-3, :], latf[:-3, :], toplot, transform=projin)
cbax = axgr.cbar_axes[cpt]
cb = cbax.colorbar(cs)
axgr[cpt].add_feature(cfeature.LAND)
axgr[cpt].add_feature(cfeature.COASTLINE)
axgr[cpt].set_extent([130, -60, -40, 40], crs=projout)
cmax = find_limits(toplot)
cs.set_clim(-cmax, cmax)
axgr[cpt].set_title('Pred.')
cb.set_label('J/m2')

plt.suptitle('L = %dcm' %l, y=0.9)
plt.savefig('oope_changes_maps_trends_l_%d' %l, bbox_inches='tight')
# -


