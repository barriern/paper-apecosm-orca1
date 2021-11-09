# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
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
import numpy as np

gear = 'PS'
species = 'SKJ'
res = 1
    
filename = 'clim_anoms_catch_gear_%s_species_%s_%dx%d.nc' %(gear, species, res, res)
filename
# -

# ## Reading Sardara file
#
# ### Reading climatology

data = xr.open_dataset(filename)
data

sarlon = data['lon'].values
sarlat = data['lat'].values

sarclim = data['clim'].where(data['clim'] > 0).mean(dim='month')
sarclim.plot()

# ### Reading anomalies

saranoms = data['anoms']
saranoms = saranoms.sel(time=slice('1997-10-01', '1997-12-31'))
saranoms = saranoms.mean(dim='time')
saranoms.plot()

# ## Reading Apecosm
#
# ### Reading climatology

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
length = const['length'].values * 100
const = const.rename({'wpred' : 'w'})
const = const.assign(w=length)
const

data = xr.open_dataset('../new-97/data/pacific_clim_OOPE.nc')
data = data.assign(w=length)
data

apeclim = data['OOPE'].mean(dim='month') * const['weight_step']
apeclim

mesh = xr.open_dataset('../new-97/data/pacific_mesh_mask.nc').isel(z=0)
mesh
lont = mesh['glamt'].values
latt = mesh['gphit'].values

# Now we process Apecosm outputs to use tricontour. First, we integrate over the size classes 25cm/70cm.

temp = apeclim.sel(w=slice(25, 70), drop=True).sum(dim='w')
temp = temp.where(temp != 0).to_masked_array()

# Now we extract the data mask array.

mask = np.ma.getmaskarray(temp)
cs = plt.imshow(mask)
plt.colorbar(cs)

# Then, we create a function that prepares the arrays for the tricontour method.

projin = ccrs.PlateCarree()
projout = ccrs.PlateCarree(central_longitude=180)


def prepare_tri(temp):

    lon1d = np.ravel(lont[~mask])
    lat1d = np.ravel(latt[~mask])
    temp1d = np.ravel(temp[~mask])
    output = projout.transform_points(projin, lon1d, lat1d)
    lonout = output[..., 0]
    latout = output[..., 1]
    
    return lonout, latout, temp1d


lonout, latout, temp1d = prepare_tri(temp)

# ## Plotting climatology

# +
plt.figure(figsize=(12, 12), facecolor='white')

ax = plt.axes(projection=projout)
step = 0.5
cl = ax.tricontour(lonout, latout, np.log10(temp1d), colors='k', linewidths=0.5, levels=np.arange(0, 3 + step, step))
plt.clabel(cl)

cs = ax.pcolormesh(sarlon, sarlat, np.log10(sarclim.values), transform=projin, shading='auto')
cb = plt.colorbar(cs, shrink=0.4)
cb.set_label('%s catch (log10)' %species)
cs.set_clim(-6, -1.5)
ax.add_feature(cfeature.LAND, zorder=2)
ax.coastlines(zorder=3)
plt.savefig('apecosm_sardara_clim_gear_%s_species_%s_%dx%d.png' %(gear, species, res, res), bbox_inches='tight')
# -
# ## Plotting anomalies


data = xr.open_dataset('../new-97/data/pacific_OOPE_anom.nc').sel(time=slice('1997-10-01', '1997-12-31'))
data = data.assign(w=length)
data

apeanom = data['OOPE'] * const['weight_step']

apeanom = apeanom.sel(w=slice(25, 70), drop=True).sum('w').mean(dim='time')
apeanom = apeanom.where(apeanom != 0)

lonout, latout, temp1d = prepare_tri(apeanom.to_masked_array())

# +
plt.figure(figsize=(12, 12), facecolor='white')

ax = plt.axes(projection=projout)
step = 50
cl = ax.tricontour(lonout, latout, temp1d, colors='k', linewidths=0.5, levels=np.arange(-300, 300 + step, step))
cl0 = ax.tricontour(lonout, latout, temp1d, colors='k', linewidths=1.2, levels=0)

cs = ax.pcolormesh(sarlon, sarlat, saranoms, transform=projin, shading='auto')
cb = plt.colorbar(cs, shrink=0.4)
#cb.set_label('%s catch (log10)' %species)
ccc = 0.01
cs.set_clim(-ccc, ccc)
ax.add_feature(cfeature.LAND, zorder=2)
ax.coastlines(zorder=3)
plt.savefig('apecosm_sardara_anoms_gear_%s_species_%s_%dx%d.png' %(gear, species, res, res), bbox_inches='tight')
# -


