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
surf = e1t * e2t * tmask
surf = np.squeeze(surf)
weights = np.sqrt(surf)
nlat, nlon = surf.shape
lon = mesh['glamt'].values
lat = mesh['gphit'].values
lon = np.squeeze(lon)
lat = np.squeeze(lat)
e3t = mesh['e3t_0'].to_masked_array() * mesh['tmask'].values
z1d = mesh['gdept_1d'].values[0]

latmax = 2
ilat, ilon = np.nonzero(np.abs(lat) < latmax)
ilat = slice(ilat.min(), ilat.max() + 1)

iz = np.nonzero(z1d <= zmax)[0]
e3t = e3t[:, iz, :, :]
e3t = e3t[:, :, ilat, :]   

surf = surf[np.newaxis, np.newaxis, ilat, :]

print(e3t.shape)  # time, depth, lat, lon
print(surf.shape)

dirin = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/JRA_CO2/'

varlist = ['thetao', 'ZOO', 'ZOO2', 'PHY2', 'O2', 'GOC', 'PAR', 'PHY', 'so']
suffix = {}
suffix['thetao'] = 'grid_T'
suffix['so'] = 'grid_T'
suffix['NCHL'] = 'add_T'
suffix['DCHL'] = 'add_T'
suffix['NO3'] = 'add_T'
suffix['PO4'] = 'add_T'
suffix['NH4'] = 'add_T'
suffix['ZOO'] = 'ptrc_T'
suffix['ZOO2'] = 'ptrc_T'
suffix['PHY'] = 'ptrc_T'
suffix['PHY2'] = 'ptrc_T'
suffix['GOC'] = 'ptrc_T'
suffix['O2'] = 'ptrc_T'
suffix['PAR'] = 'diad_T'

varlist = ['thetao', 'so', 'ZOO', 'ZOO2', 'PHY', 'PHY2', 'O2', 'GOC', 'PAR']
    
for varname in varlist:

    print('Processing ++++++++++++++++++++++++ ', varname)

    datatot = xr.open_mfdataset("%s/*%s.nc" %(dirin, suffix[varname]), combine='by_coords')
    datatot = datatot.isel(olevel=iz, y=ilat)
    year = datatot['time_counter'].values

    data = datatot[varname]

    data = data.to_masked_array()
    data = np.sum(data * e3t * surf, axis=(1, 2)) / np.sum(e3t * surf, axis=(1, 2))
    #dataout = detrend(year, data)
    
    dataout = data
    print(dataout.shape)
    
    dsout = xr.Dataset({varname:(['year', 'x'], dataout), 
                        'year':(['year'], year)})
    dsout.attrs['creation_date'] = str(date.today())
    dsout.attrs['script'] = os.path.realpath(__file__)
    dsout.attrs['description'] = 'Detrended 0-%d anomalies of %s, |lat| < %f' %(zmax, varname, latmax)
    dsout.to_netcdf('%s/hovmoller_data_0-%d_%s.nc' %(dirout, zmax, varname))
