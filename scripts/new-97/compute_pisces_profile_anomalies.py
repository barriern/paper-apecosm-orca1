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
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from dask.diagnostics import ProgressBar

latmax = 1
grid = 'speed_U'
varname = 'uo'
cell = 'u'
dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/processed_pacific'
# -

# ## Reading mesh mask

mesh = xr.open_dataset('%s/pacific_mesh_mask.nc' %dirin)
mesh = mesh.rename({'z': 'olevel'})
mesh

lat = mesh['gphif']
lat

mesh = mesh.where(abs(lat) <= latmax)

lon0 = mesh['glamt'].mean(dim='y').isel(olevel=0)
lon0 = (lon0 + 360) % 360
lon0

maskvarname = cell + 'mask'
e3varname = 'e3' + cell + '_0'
e2varname = 'e2' + cell
e1varname = 'e1' + cell
volume = mesh[maskvarname] * mesh[e3varname] * mesh[e1varname] * mesh[e2varname]
volume.name = 'volume'
volume

# ## Loading the climatology

clim = xr.open_dataset('%s/pacific_clim_%s.nc' %(dirin, grid))
clim

varclim = clim[varname]
varclim

# ## Loading the field

data = xr.open_dataset('%s/pacific_nino97_%s.nc' %(dirin, grid))
data

ntime = data.dims['time_counter']
ntime

var = data[varname]
var

var = var.chunk({'time_counter': 12})
var

# ## Computing the anomalies

anom = var.groupby('time_counter.month') - varclim
anom

# ## Computing the mean anomalies profile

anomprofile = (anom * volume).sum(dim=['y']) / (volume.sum(dim=['y']))
anomprofile

anomprofile['x'] = lon0.values
anomprofile

# ## Saving the anomalies

outfile = '%s/pacific_profile_anoms_%s.nc' %(dirin, varname)
outfile

delayed = anomprofile.to_netcdf(outfile, compute=False)

with ProgressBar():
    delayed.compute()


