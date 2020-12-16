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
zmax = 200

# Load the mesh mask
mesh = xr.open_dataset("/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/corrected_mesh_mask_eORCA1_v2.2.nc")
tmask = np.squeeze(mesh['tmask'].values[0, 0])
e1t = mesh['e1t'].values
e2t = mesh['e2t'].values
surf = e1t * e2t
surf = np.squeeze(surf)
weights = np.sqrt(surf)
nlat, nlon = surf.shape
lon = mesh['glamt'].values
lat = mesh['gphit'].values
lon = np.squeeze(lon)
lat = np.squeeze(lat)
e3t = mesh['e3t_0'].values * mesh['tmask'].values
z1d = mesh['gdept_1d'].values[0]

iz = np.nonzero(z1d <= zmax)[0]
e3t = e3t[:, iz, :, :]

dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/corr_mask/output/yearly/pisces'

varlist = ['thetao', 'ZOO', 'ZOO2', 'PHY2', 'O2']
varlist = ['thetao', 'ZOO', 'ZOO2', 'PHY2', 'O2', 'GOC', 'PAR', 'PHY', 'so']
varlist = ['CHL']
varlist = ['uo', 'vo']
    
datatot = xr.open_mfdataset("%s/yearly_mean_*.nc" %(dirin), combine='by_coords')
datatot = datatot.isel(olevel=iz)
year = datatot['time_counter'].values

for varname in varlist:

    print('Processing %s' %varname)
    if varname == 'CHL':
        data = datatot['NCHL'] + datatot['DCHL']
    elif varname == 'PLK':
        data = datatot['ZOO'] + datatot['ZOO2'] + datatot['PHY'] + datatot['PHY2'] + datatot['GOC']
    else:
        data = datatot[varname]

    data = data.to_masked_array()
    data = np.sum(data * e3t, axis=1) / np.sum(e3t, axis=1)
    dataout = detrend(year, data)

    dsout = xr.Dataset({varname:(['year', 'y', 'x'], dataout), 
                        'year':(['year'], year)})
    dsout.attrs['creation_date'] = str(date.today())
    dsout.attrs['script'] = os.path.realpath(__file__)
    dsout.attrs['description'] = 'Detrended 0-%d anomalies of %s' %(zmax, varname)
    dsout.to_netcdf('%s/detrended_0-%d_yearly_%s.nc' %(dirout, zmax, varname))
