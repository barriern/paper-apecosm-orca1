# # Computation of simulated ONI index
#
# The aim of this script is to compute the simulated TPI index
#
# ## Import libraries

import numpy as np
from datetime import date
import os.path
import numpy as np
from eofs.standard import Eof
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import stats
import xarray as xr
from glob import glob
import apecosm.ts as ts

# ## Defining paths
#
# First, we define the input directory and the file to process:

dirin = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/JRA_CO2/'
filelist = glob('%s/*grid_T*nc' %dirin)
filelist.sort()
filelist

# Then, we define the output directory:

dirout = '../data'

# ## Reading the mesh mask

mesh = xr.open_dataset("../../../data/mesh_mask_eORCA1_v2.2.nc").isel(t=0)
tmask = np.squeeze(mesh['tmask'].values[0])
e1t = mesh['e1t'].values
e2t = mesh['e2t'].values
surf = e1t * e2t
surf = np.squeeze(surf)
nlat, nlon = surf.shape
lon = mesh['glamt'].values
lat = mesh['gphit'].values
lon = np.squeeze(lon)
lat = np.squeeze(lat)
lonf = mesh['glamf'].values
latf = mesh['gphif'].values

# Now, the ONI domain (Nino 3.4) domain is extracted, as indicated in https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni

latmin = -5
latmax = 5
lonmax = -120
lonmin = -170

test = (lat <= latmax) & (lat >= latmin)
test = test & (lon<=lonmax) & (lon>=lonmin)
test = test & (tmask == 1)
xrsurf = mesh['e1t'] * mesh['e2t'] * mesh['tmaskutil']
xrsurf = xrsurf.where(test)
xrsurf.plot()

ilat, ilon = np.nonzero(test == True)

# Now, we make a plot to see if the right domain is extracted:

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

output = tmask.copy()
output[ilat, ilon] = 2

plt.figure()
ax = plt.axes(projection=proj)
cs = ax.pcolormesh(lonf, latf, output[1:, 1:], transform=proj2)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
plt.colorbar(cs)
plt.show()


# ## Processing SST files
#
# Now, the NEMO SST files are processed. First, the climatology between January 1971 and December 2000 is extracted.

# ### Computation of climatology

data = xr.open_mfdataset(filelist).isel(olevel=0)
data = data['thetao']
years = data['time_counter.year'].values
months = data['time_counter.month'].values
date = years * 100 + months
data

clim = data.sel(time_counter=slice('1971-01-01', '2000-12-31'))
clim = clim.groupby('time_counter.month').mean(dim='time_counter')
clim

# ### Computation of anomalies

anom = data.groupby('time_counter.month') - clim
anom

# ### Computation of ONI index
#
# Now the weighted mean of the SST anomalies is computed, with the weights provided by the cell surface.

output = (anom * xrsurf).sum(dim=['x', 'y']) / xrsurf.sum(dim=['x', 'y'])
output.name = 'enso'
output

output.to_netcdf('../data/simulated_enso_index.nc')
