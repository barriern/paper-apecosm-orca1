# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt

zmax = 200
varname = 'GOC'
latmax = 1
# -

# ## Reading mesh mask

mesh = xr.open_dataset('data/pacific_mesh_mask.nc')
mesh

lat = mesh['gphit']
lat

mesh = mesh.where(abs(lat) <= latmax)

lon = mesh['glamt'].mean(dim='y').isel(z=0)
lon

lon = (lon + 360) % 360
lon

surf = mesh['e1t'] * mesh['e2t'] * mesh['tmask']
surf = surf.isel(z=0)
surf

# ## Loading the climatology

clim = xr.open_dataset('data/pacific_0-%d_%s_clim.nc' %(zmax, varname))
clim

varclim = clim[varname]
varclim

# ## Loading the field

data = xr.open_dataset('data/pacific_0-%d_%s.nc' %(zmax, varname))
data

var = data[varname]
var

# ## Computing the anomalies

anom = var.groupby('time_counter.month') - varclim
anom

# ## Computing the meridional anomalies

anomout = (anom * surf).sum(dim='y') / surf.sum(dim='y')
anomout

anomout['x'] = lon

plt.figure(facecolor='white')
anomout.plot(robust=True)
plt.gca().set_title('0-%d %s' %(zmax, varname))
figname = 'figs/pacific_0-%d_anom_hov_%s.png' %(zmax, varname)
plt.savefig(figname)


