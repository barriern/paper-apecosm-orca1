import xarray as xr
import numpy as np
import os.path 
from datetime import datetime
from glob import glob

prefix = 'final-runs'

dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/%s/output/subclass' %(prefix)
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'
dirmesh = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST'

latmax = 2

mesh = xr.open_dataset('%s/corrected_mesh_mask_eORCA1_v2.2.nc' %dirmesh)
lon = mesh['glamt'].values[0]
lat = mesh['gphit'].values[0]

ilat, ilon = np.nonzero(np.abs(lat) < latmax)
jmin = ilat.min()
jmax = ilat.max() + 1

mesh = mesh.isel(y=slice(jmin, jmax))
lon = mesh['glamt'].values[0]
lat = mesh['gphit'].values[0]
tmask = mesh['tmask'].values[0, 0]
e1t = mesh['e1t'].values[0, 0]
e2t = mesh['e2t'].values[0, 0]
surf = e1t * e2t * tmask
surf = surf[np.newaxis, :, :, np.newaxis, np.newaxis]  

lon0 = np.mean(lon, axis=0)

for varname in ['diff', 'mdiff_trend', 'zdiff_trend', 'starvation', 'u_active', 'u_passive', 'v_active', 'v_passive', 'madv_trend', 'zadv_trend']:
#for varname in ['OOPE']:

    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@", varname)
   
    pattern = '%s/*%s*nc' %(dirin, varname)

    data = xr.open_mfdataset(pattern, combine='by_coords')
    data = data.isel(y=slice(jmin, jmax))
    data = data[varname] # time, lat, lon, com, size (732, 11, 150, 3, 3)

    print(surf.shape, data.shape)

    dataout = (data * surf).sum(dim='y') / (np.sum(surf, axis=1))
    dataout.coords['x'] = lon0

    fileout = '%s/%s_%s_meridional_mean.nc' %(dirout, prefix, varname)
    print(fileout)
    dataout.to_netcdf(fileout)
    dataout.attrs['file'] = os.path.realpath(__file__)
    dataout.attrs['date'] = str(datetime.today())
