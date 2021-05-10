import xarray as xr
import numpy as np
import os.path 
from datetime import datetime
from glob import glob
import cftime
from datetime import datetime

units = 'seconds since 1950-01-01 00:00:00'
calendar='noleap'

latmax = 2

dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data/'

dirmesh = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST'
mesh = xr.open_dataset('%s/corrected_mesh_mask_eORCA1_v2.2.nc' %dirmesh)
lon = mesh['glamt'].values[0]
lat = mesh['gphit'].values[0]
e1t = mesh['e1t'].values[0]
e2t = mesh['e2t'].values[0]
tmask = mesh['tmask'].values[0, 0]
surf = e1t * e2t * tmask  # lat, lon

ilat, ilon = np.nonzero(np.abs(lat) < latmax)
jmin = ilat.min()
jmax = ilat.max() + 1
ilat = slice(jmin, jmax)

surf = surf[ilat, :]
lon0 = np.mean(lon[ilat, :], axis=0)
surf = surf[np.newaxis, np.newaxis, :, :, np.newaxis, np.newaxis, np.newaxis]   # time, dn, y, x, depth, com, size

prefix = 'final-runs'
dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/%s/output' %prefix

years = np.arange(1958, 2019)

months = np.arange(12) + 1

filelist = np.sort(glob('%s/*FORAGE*nc' %dirin))

for f in filelist:

    bname = os.path.basename(f)

    data = xr.open_mfdataset(f, combine='by_coords')
    data = data.isel(y=ilat)
    forage = data['FORAGE']  # time, dn, y, x, depth, com, size
    forage = forage.where(forage != -999)  # mask filled values
    forage = forage * surf
    forage = forage.sum('y') / np.sum(surf, axis=2)
    #forage.coords['time'] = time
    #forage.coords['time'].attrs['units'] = units
    #forage.coords['time'].attrs['calendar'] = calendar

    fout = '%s/equatorial_%s' %(dirout, bname)

    forage.attrs['file'] = os.path.realpath(__file__)
    
    forage.to_netcdf(fout, unlimited_dims='time')
