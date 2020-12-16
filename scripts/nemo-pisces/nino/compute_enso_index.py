import numpy as np
from datetime import date
import os.path
import numpy as np
from eofs.standard import Eof
import matplotlib.pyplot as plt
import cartopy.crs as crs
from scipy import stats
import xarray as xr
from glob import glob

dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'
dirin = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/JRA_CO2/'

# Load the mesh mask
mesh = xr.open_dataset("/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/mesh_mask_eORCA1_v2.2.nc")
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

surf = surf * tmask
surf = surf[np.newaxis, :, :] 

latmin = -5
latmax = 5
lonmax = -120
lonmin = -170

test = (lat<=latmax) & (lat>=latmin)
test = test & (lon<=lonmax) & (lon>=lonmin)
test = test & (tmask == 1)

ilat, ilon = np.nonzero(test == True)
surf = surf[:, ilat, ilon]

filelist = np.sort(glob('%s/*grid_T*nc' %dirin))

for f in filelist:

    print('processing %s' %f)

    data = xr.open_dataset('%s' %(f))
    data = data.isel(olevel=0)
    time = data['time_counter']
    data = data['thetao'].values

    ntime, nlat, nlon = data.shape

    temp = data[:, ilat, ilon] * surf
    temp = np.sum(temp, axis=1) / np.sum(surf, axis=1)

    if f == filelist[0]:
        output = temp
        timeout = time
    else:
        output = np.concatenate((output, temp), axis=0)
        timeout = np.concatenate((timeout, time), axis=0)

dsout = xr.Dataset({'enso':(['time'], output)})
dsout['time'] = (['time'], timeout)
dsout.attrs['creation_date'] = str(date.today())
dsout.attrs['script'] = os.path.realpath(__file__)
dsout.attrs['description'] = 'EOF for densities by class'
dsout.to_netcdf('%s/simulated_enso_index.nc' %(dirout))
