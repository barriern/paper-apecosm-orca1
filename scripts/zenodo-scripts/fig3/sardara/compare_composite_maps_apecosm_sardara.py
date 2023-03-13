# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python [conda env:nbarrier]
#     language: python
#     name: conda-env-nbarrier-py
# ---

# +
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np

ilat = slice(None, -3)
# -

mesh = xr.open_dataset('../data/static/pacific_mesh_mask.nc').isel(z=0, y=ilat)
mesh
lonf = mesh['glamf'].values
latf = mesh['gphif'].values
lonf

const = xr.open_dataset('../data/static/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100
const

ape_nina = xr.open_dataset('compo_apecosm_nina.nc').isel(y=ilat)
ape_nina = ape_nina.rename({'w': 'l'})
ape_nina['l'] = const['l']
ape_nina = (ape_nina * const['weight_step']).sel(l=slice(30, 70)).sum(dim='l')
ape_nina

ape_nino = xr.open_dataset('compo_apecosm_nino.nc').isel(y=ilat)
ape_nino = ape_nino.rename({'w': 'l'})
ape_nino['l'] = const['l']
ape_nino = (ape_nino * const['weight_step']).sel(l=slice(30, 70)).sum(dim='l')
ape_nino

ape_diff = ape_nino - ape_nina

sar_nina = xr.open_dataset('compo_sardara_nina.nc')
sar_nina

sar_nino = xr.open_dataset('compo_sardara_nino.nc')
sar_nino

sar_diff = sar_nino - sar_nina
sar_diff = sar_diff.where(sar_diff != 0)

rt = 6371 * 1e3
dy = rt * np.deg2rad(1)
dx = rt * np.cos(np.deg2rad(sar_diff['lat'])) * np.deg2rad(1)
surface = dx * dy
surface
sar_diff = sar_diff / surface
sar_diff

ape_diff = ape_nino - ape_nina
ape_diff = ape_diff.where(ape_diff != 0)
ape_diff /=  (4e6 * 1000)

# +
plt.figure(figsize=(12,8), facecolor='white')
plt.rcParams['font.size'] = 15
plt.rcParams['image.cmap'] = 'RdBu_r'


projout = ccrs.PlateCarree(central_longitude=180)
projin = ccrs.PlateCarree()

ax = plt.subplot(2, 1, 1, projection=projout)
cs = plt.pcolormesh(sar_diff['lon'].values, sar_diff['lat'].values, sar_diff['catch'].values, shading='auto', transform=projin)
cb = plt.colorbar(cs)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
ccc = 3e-8
cs.set_clim(-ccc, ccc)
cb.set_label('Catch (Ton/m2)')
plt.title('Sardara')
ax.set_extent([130, -60 + 360, -40, 40], crs=projin)



ax = plt.subplot(2, 1, 2, projection=projout)
cs = plt.pcolormesh(lonf, latf, ape_diff['OOPE'].values[1:, 1:], transform=projin)
ccc = 4e-8
cb = plt.colorbar(cs)
cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
cb.set_label('Biomass (Ton/m2)')
plt.title('Apecosm')
ax.set_extent([130, -60 + 360, -40, 40], crs=projin)
plt.savefig('composites_nino-nina', bbox_inches='tight')
# -

ape_diff.to_netcdf('map_to_plot_ape.nc')
sar_diff.to_netcdf('map_to_plot_sar.nc')


