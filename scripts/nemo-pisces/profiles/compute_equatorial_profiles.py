''' Extracts and averages Pisces data over the 0-200m layer '''

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
latmax = 2

# Load the mesh mask
mesh = xr.open_dataset("/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/corrected_mesh_mask_eORCA1_v2.2.nc")
tmask = np.squeeze(mesh['tmask'].values[0])  # z, y, x
e1t = mesh['e1t'].values
e2t = mesh['e2t'].values
surf = e1t * e2t   
surf = np.squeeze(surf) * tmask  # z, y, z
nz, nlat, nlon = surf.shape
lon = mesh['glamt'].values
lat = mesh['gphit'].values
lon = np.squeeze(lon)
lat = np.squeeze(lat)

ilat, ilon = np.nonzero(np.abs(lat) <= latmax)
jmin, jmax = ilat.min(), ilat.max() + 1

ilat = slice(jmin, jmax)

lon0 = np.mean(lon[ilat, :], axis=0)
surf = surf[np.newaxis, :, ilat, :]  # 1, z, lat, lon

surf_int = np.sum(surf, axis=2)

dirin = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/JRA_CO2/'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data' 

filelist = np.sort(glob('%s/*ptrc_T*nc' %dirin))
#filelist = np.sort(glob('%s/nico*1958*grid_W*nc' %dirin))
print(filelist)

for f in filelist:

    if('ssh' in f):
        continue
    if('grid_W' in f):
        continue
    
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ File %s' %f)

    fout = '%s/equatorial_%s' %(dirout, os.path.basename(f))

    data = xr.open_dataset(f)
    data = data.isel(y=ilat)

    var = data.variables
    varout = [v for v in var if ((data[v].ndim == 4) & ('e3' not in v))]

    dsout = xr.Dataset()
    for v in varout:
        print('Variable ', v)
        print(data[v].shape)  # 12, 75, 11, 362
        print(surf.shape)  # 1, 75, 11, 362
        print((data[v] * surf).shape)  #  12,75, 11, 362
        print(surf_int.shape)  # 1, 75, 1, 362
        # num = (12, 75, 362)
        temp = (data[v] * surf).sum(dim='y') / surf_int
        dsout[v] = temp
    
    dsout.attrs['file'] = os.path.realpath(__file__)
    dsout.attrs['des2'] = 'Average between -%d and %d lat' %(latmax, latmax)
    dsout.to_netcdf(fout, unlimited_dims='time_counter')
