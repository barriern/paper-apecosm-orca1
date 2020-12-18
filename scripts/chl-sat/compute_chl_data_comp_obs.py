import numpy as np
from datetime import date
import os.path
import numpy as np
from eofs.standard import Eof
import matplotlib.pyplot as plt
import cartopy.crs as crs
from scipy import stats
import xarray as xr
from remove_trend_yearly import detrend

dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

# Load the mesh mask
mesh = xr.open_dataset("/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/corrected_mesh_mask_eORCA1_v2.2.nc")
mesh = mesh.isel(z=0)
tmask = np.squeeze(mesh['tmask'].values[0])
e1t = mesh['e1t'].values
e2t = mesh['e2t'].values
surf = e1t * e2t * tmask
surf = np.squeeze(surf)
weights = np.sqrt(surf)
nlat, nlon = surf.shape
lon = mesh['glamt'].values
lat = mesh['gphit'].values
lon = np.squeeze(lon)
lat = np.squeeze(lat)

latmax = 2
ilat, ilon = np.nonzero(np.abs(lat) < latmax)
ilat = slice(ilat.min(), ilat.max() + 1)

surf = surf[np.newaxis, ilat, :]
lon0 = np.mean(lon[ilat, :], axis=0)



dirin = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/JRA_CO2/'
pattern = "%s/*add_T.nc" %dirin

data = xr.open_mfdataset(pattern, combine='by_coords')
data = data.isel(y=ilat, olevel=0)
chl = data['NCHL'] + data['DCHL']

print(chl.shape)
print(surf.shape)
output = (chl * surf).sum(dim='y') / np.sum(surf, axis=1)
output = output.rename('chl')
output['x'] = lon0
output.to_netcdf('%s/simulated_surface_chl_hov.nc' %dirout)
