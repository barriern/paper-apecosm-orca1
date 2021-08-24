# # Computation of simulated TPI index
#
# Adapted from https://psl.noaa.gov/data/timeseries/IPOTPI/tpi.create.ncl
#
# ## Importing modules

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

dirin = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/JRA_CO2'
os.listdir(dirin)
filelist = glob('%s/surface*nc' %dirin)
filelist.sort()
filelist

# Then, we define the output directory:

dirout = os.getenv('DATAWORK')
dirout = os.path.join(dirout, 'apecosm', 'apecosm_orca1', 'diags', 'final_diags')
dirout


# ## Extraction of the domain
#
# Creation of a function to extract a given domain. It checks latitudes boundaries and `tmask`.

def extract_domain(latmin, latmax, lonmin, lonmax):
    test = ((lat >= latmin) & (lat <= latmax))
    test2 = (lon >= lonmin) | (lon <= lonmax)
    test = (test & test2) * tmask
    iok1 = np.nonzero(test)
    return iok1


# Note that the management of longitudes is not done in the same way as for latitude, since the focus is on the Pacific basin, where date change line is located. Therefore, $lonmin>0$ while $lonmax<0$.

# ## Extraction of the mesh mask

mesh = xr.open_dataset("../../../data/mesh_mask_eORCA1_v2.2.nc")
tmask = np.squeeze(mesh['tmask'].values[0, 0])
e1t = mesh['e1t'].values
e2t = mesh['e2t'].values
surf = e1t * e2t
surf = np.squeeze(surf)
nlat, nlon = surf.shape
lon = mesh['glamt'].values
lat = mesh['gphit'].values
lon = np.squeeze(lon)
lat = np.squeeze(lat)
lonf = mesh['glamf'].values[0]
latf = mesh['gphif'].values[0]

# ## Definition of the mask used for the computation

# The variable that is used to define the domain mask is based on `tmask`:
# - $0$: land
# - $1$: water but out of domain
# - $100$: first domain
# - $200$: second domain
# - $300$: third domain
#
# First, a deep copy of the original tmask variable is made.

output = tmask.copy()

# Extraction of the first domain:

latmin = 25
latmax= 45
lonmin = 140
lonmax = -145
iok1 = extract_domain(latmin, latmax, lonmin, lonmax)
output[iok1] = 100

# Extraction of the second domain:

latmin = -10
latmax= 10
lonmin = 170
lonmax = -90
iok2 = extract_domain(latmin, latmax, lonmin, lonmax)
output[iok2] = 200

# Extraction of the third domain:

latmin = -50
latmax= -15
lonmin = 150
lonmax = -160
iok3 = extract_domain(latmin, latmax, lonmin, lonmax)
output[iok3] = 300

# Now, the domain definition is verified, in order to insure that the proper means are computed.

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

plt.figure()
ax = plt.gca(projection=proj)
cs = ax.pcolormesh(lonf, latf, output[1:, 1:], transform=proj2)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
plt.colorbar(cs)
plt.savefig('boxes.png')
plt.show()

# ## Processing of the SST values
#
# Now that the domain are defined, the SST data are processed.
#
# ### Computation of climatology
#
# Climatology is computed between 1971 and 2000. 

datestart = 197101
dateend = 200012

data = xr.open_mfdataset(filelist).isel(olevel=0)
years = data['time_counter.year'].values
months = data['time_counter.month'].values
date = years * 100 + months
data

iclim = np.nonzero((date >= datestart) & (date <= dateend))[0]
iclim

sst = np.squeeze(data['thetao'].values)
ntime, nlat, nlon = sst.shape

clim, anom = ts.get_monthly_clim(sst[iclim, :, :])

# ### Computation of SST anomalies

nyears = ntime // 12
index = np.arange(12)

anom = np.zeros(sst.shape)
for y in range(nyears):
    anom[index] = sst[index] - clim
    index += 12

# ### Computation of TPI index

# First, the surface ($e1t \times e2t$) variable is modified to include a spurious time dimension.

surf = surf[np.newaxis, :]
surf.shape

# First, computation of mean SST anomalies over the first domain

ilat = iok1[0]
ilon = iok1[1]
ts1 = np.sum(surf[:, ilat, ilon] * anom[:, ilat, ilon], axis=-1) / np.sum(surf[:, ilat, ilon], axis=-1)

# Computation of mean SST anomalies over the second domain

ilat = iok2[0]
ilon = iok2[1]
ts2 = np.sum(surf[:, ilat, ilon] * anom[:, ilat, ilon], axis=-1) / np.sum(surf[:, ilat, ilon], axis=-1)

# Second, computation of mean SST anomalies over the third domain

ilat = iok3[0]
ilon = iok3[1]
ts3 = np.sum(surf[:, ilat, ilon] * anom[:, ilat, ilon], axis=-1) / np.sum(surf[:, ilat, ilon], axis=-1)

# Finally, computation of TPI index.

tpi = ts2 - 0.5 * (ts1 + ts3)

# ## Saving outputs

# + language="javascript"
# IPython.notebook.kernel.execute('nb_name = "' + IPython.notebook.notebook_name + '"')
# -

__file__ = os.path.join(os.getcwd(), nb_name)
__file__

dsout = xr.Dataset()
dsout['ts1'] = (['time'], ts1)
dsout['ts2'] = (['time'], ts2)
dsout['ts3'] = (['time'], ts3)
dsout['tpi'] = (['time'], tpi)
dsout['time'] = date
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('model_tpi_index.nc')
