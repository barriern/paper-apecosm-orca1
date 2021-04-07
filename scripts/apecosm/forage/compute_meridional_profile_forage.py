import xarray as xr
import numpy as np
import os.path 
from datetime import datetime
from glob import glob
import cftime
from datetime import datetime

units = 'seconds since 1950-01-01 00:00:00'
calendar='noleap'

lonmax = -150

dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data/'

dirmesh = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST'
mesh = xr.open_dataset('%s/corrected_mesh_mask_eORCA1_v2.2.nc' %dirmesh)
lon = mesh['glamt'].values[0]
lat = mesh['gphit'].values[0]
e1t = mesh['e1t'].values[0]
e2t = mesh['e2t'].values[0]
surf = e1t * e2t   # lat, lon

# extraction of the proper longitude indexes
ilat, ilon = np.nonzero(np.abs(lat) <= 20)
jmin = ilat.min()
jmax = ilat.max() + 1
ilat = slice(jmin, jmax)
lon0 = np.mean(lon[ilat, :], axis=0)


ilon = np.nonzero((lon0 <= lonmax + 2) & (lon0 >= lonmax - 2))[0]
imin, imax = ilon.min(), ilon.max() + 1
ilon = slice(imin, imax)

lat0 = np.mean(lat[:, ilon], axis=-1)

surf = surf[:, ilon]  # lat, lon
surf = surf[np.newaxis, np.newaxis, :, :, np.newaxis, np.newaxis, np.newaxis]  # time, dn, y, x, depth, w

prefix = 'final-runs'
dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/%s/output' %prefix

years = np.arange(1958, 2019)

months = np.arange(12) + 1

print(surf.shape)

for y in years:

    date = [datetime(y, m, 15) for m in months]
    date = np.array(date)
    time = cftime.date2num(date, units, calendar)

    pattern = '%s/*FORAGE_Y%.04dD030.nc' %(dirin, y)
    print(pattern)
    filelist = np.sort(glob(pattern))
    data = xr.open_mfdataset(filelist, combine='by_coords')
    data = data.isel(x=ilon)
    forage = data['FORAGE']  # time, dn, y, x, depth, size
    #forage = forage.where(forage != -999)  # mask filled values
    print(forage.shape)
    forage = forage * surf
    
    forage = forage.sum('x') / np.sum(surf, axis=3)
    forage.coords['time'] = time
    forage.coords['time'].attrs['units'] = units
    forage.coords['time'].attrs['calendar'] = calendar
    forage['y'] = lat0
   
    fout = 'meridional_%.f_forage_year_%.4d.nc' %(lonmax, y)
    fout = '%s/%s' %(dirout, fout)
    forage.attrs['file'] = os.path.realpath(__file__)
    forage.to_netcdf(fout, unlimited_dims='time')
