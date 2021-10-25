# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar

plt.rcParams['text.usetex'] = False

zmax = 50
grid = 'ptrc_T'
e3var = 'e3t_0'
maskvar = 'tmask'
varname = 'ZOO2'
# -

# ## Extracting the mesh mask

mesh = xr.open_dataset('data/pacific_mesh_mask.nc')
mesh

z1d = mesh['gdept_1d']
z1d.shape

mesh = mesh.where(z1d <= zmax)
mesh['gdept_1d'].isel(x=45, y=45).values

mesh = mesh.rename({'z': 'olevel'})
mesh

e3t = mesh[e3var]
e3t

mask = mesh[maskvar]
mask

# ## Extracting the variable

data = xr.open_dataset('data/pacific_nino97_%s.nc' %grid)
data

var = data[varname]
var

varout = (var * e3t * mask).sum(dim='olevel') / (e3t * mask).sum(dim='olevel')
varout.name = varname
varout

outfile = 'data/pacific_0-%s_%s.nc' %(zmax, varname)
outfile

delayed = varout.to_netcdf(outfile, compute=False)
with ProgressBar():
    delayed.compute()

# ## Processing climatology

clim = xr.open_dataset('data/pacific_clim_%s.nc' %grid)
clim

varclim = clim[varname]
varclim

varclimout = (varclim * e3t * mask).sum(dim='olevel') / (e3t * mask).sum(dim='olevel')
varclimout.name = varname
varclimout

outfile = 'data/pacific_0-%s_%s_clim.nc' %(zmax, varname)
outfile

delayed = varclimout.to_netcdf(outfile, compute=False)
with ProgressBar():
    delayed.compute()


