''' Scripts that extracts and detrend a zonal profil from yearly Pisces data. '''
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

# output directory
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'
zmax = 1000

# Load the mesh mask
mesh = xr.open_dataset("/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/corrected_mesh_mask_eORCA1_v2.2.nc")
tmask = np.squeeze(mesh['tmask'].values[0, 0])
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
e3t = mesh['e3t_0'].values * mesh['tmask'].values
z1d = mesh['gdept_1d'].values[0]

# extracts the latitude index between -2/+2
latmax = 2
ilat, ilon = np.nonzero(np.abs(lat) < latmax)
jmin, jmax = ilat.min(), ilat.max()

# extracts the surface between -2/+2 and averages longitude between these lat
surf = surf[np.newaxis, np.newaxis, jmin:jmax+1, :]
lon0 = lon[jmin:jmax+1, :]
lon0 = np.mean(lon0, axis=0)

# extract e3t between them
e3t = e3t[:, :, jmin:jmax + 1, :]

dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/corr_mask/output/yearly/pisces'

varlist = ['thetao', 'ZOO', 'ZOO2', 'PHY2', 'O2', 'GOC', 'PAR', 'PHY', 'uocetr_eff', 'vocetr_eff', 'so', 'CHL', 'PLK']
varlist = ['uo', 'vo']
    
datatot = xr.open_mfdataset("%s/yearly_mean_*.nc" %(dirin), combine='by_coords')
year = datatot['time_counter'].values

for varname in varlist:

    print('Processing %s' %varname)
    if varname == 'CHL':
        data = datatot['NCHL'] + datatot['DCHL']
    elif varname == 'PLK':
        data = datatot['ZOO'] + datatot['ZOO2'] + datatot['PHY'] + datatot['PHY2'] + datatot['GOC']
    else:
        data = datatot[varname]

    # extracts the data
    data = data.to_masked_array()
    data = data[:, :, jmin:jmax+1, :]  # time, depth, lat, lon
    
    # computes volume weighted mean along latitude
    data = np.sum(data * surf * e3t, axis=2) / np.sum(surf * e3t, axis=2)
    
    # detrend data
    dataout = detrend(year, data)

    dsout = xr.Dataset({varname:(['year', 'depth', 'lon'], dataout), 
                        'lon':(['lon'], lon0),
                        'depth':(['depth'], z1d),
                        'year':(['year'], year)})
    dsout.attrs['creation_date'] = str(date.today())
    dsout.attrs['script'] = os.path.realpath(__file__)
    dsout.attrs['description'] = 'Detrended profile anomalies of %s, |lat| < %f' %(varname, latmax)
    dsout.to_netcdf('%s/detrended_profile_yearly_%s.nc' %(dirout, varname))
