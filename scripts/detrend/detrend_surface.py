''' Scripts that extracts Pisces yearly data at the surface, and detrend it '''

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
tmask = np.squeeze(mesh['tmask'].values[0, 0])
lon = mesh['glamt'].values
lat = mesh['gphit'].values
lon = np.squeeze(lon)
lat = np.squeeze(lat)

dirin =  '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/yearly_pisces'

# List of variables to extract
varlist = ['thetao', 'ZOO', 'ZOO2', 'PHY2', 'O2']
varlist = ['PAR', 'PHY', 'so']
varlist = ['thetao']

for varname in varlist:

    print('Processing %s' %varname)

    data = xr.open_dataset("%s/yearly_mean_%s.nc" %(dirin, varname))
    year = data['time_counter'].values
    data = data.isel(olevel=0)
    data = data[varname].to_masked_array()

    dataout = detrend(year, data)

    dsout = xr.Dataset({varname:(['year', 'y', 'x'], dataout), 
                        'year':(['year'], year)})
    dsout.attrs['creation_date'] = str(date.today())
    dsout.attrs['script'] = os.path.realpath(__file__)
    dsout.attrs['description'] = 'Detrended surface anomalies of %s' %varname
    dsout.to_netcdf('%s/detrended_surf_yearly_%s.nc' %(dirout, varname))
