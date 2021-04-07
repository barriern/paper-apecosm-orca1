import xarray as xr
import numpy as np
import os.path 
from datetime import datetime
from glob import glob

prefix = 'final-runs'

dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/%s/output/subclass' %(prefix)
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'
dirmesh = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST'

latmax = 20
lonmax = -150

mesh = xr.open_dataset('%s/corrected_mesh_mask_eORCA1_v2.2.nc' %dirmesh)
lon = mesh['glamt'].values[0]
lat = mesh['gphit'].values[0]

ilat, ilon = np.nonzero(np.abs(lat) <= latmax)
jmin = ilat.min()
jmax = ilat.max() + 1

lon0 = lon[slice(jmin, jmax), :].mean(axis=0)
ilon = np.nonzero((lon0 <= lonmax + 2) & (lon0 >= lonmax - 2))[0]
imin = ilon.min()
imax = ilon.max() + 1

mesh = mesh.isel(x=slice(imin, imax))
lon = mesh['glamt'].values[0]
lat = mesh['gphit'].values[0]
tmask = mesh['tmask'].values[0, 0]
e1t = mesh['e1t'].values[0, 0]
e2t = mesh['e2t'].values[0, 0]
surf = e1t * e2t * tmask
surf = surf[np.newaxis, :, :, np.newaxis, np.newaxis]   # time, lat, lon, c, w 

lat0 = np.mean(lat, axis=1) 

for varname in ['OOPE', 'diff', 'mdiff_trend', 'zdiff_trend', 'starvation', 'u_active', 'u_passive', 'v_active', 'v_passive', 'madv_trend', 'zadv_trend']:
#for varname in ['gamma1', 'mort_day']:

    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@", varname)
   
    pattern = '%s/*%s*nc' %(dirin, varname)
    print(pattern)

    try:
        data = xr.open_mfdataset(pattern, combine='by_coords')
        data = data.isel(x=slice(imin, imax))
        data = data[varname] # time, lat, lon, com, size (732, 11, 150, 3, 3)

        print(surf.shape, data.shape)

        dataout = (data * surf).sum(dim='x') / (np.sum(surf, axis=2))
        dataout.coords['y'] = lat0

        fileout = '%s/zonal_mean_%.f_%s_%s.nc' %(dirout, lonmax, prefix, varname)
        print(fileout)
        dataout.to_netcdf(fileout)
        dataout.attrs['file'] = os.path.realpath(__file__)
        dataout.attrs['date'] = str(datetime.today())
    except:
        pass
