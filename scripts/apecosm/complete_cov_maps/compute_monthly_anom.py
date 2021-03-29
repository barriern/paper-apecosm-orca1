import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime
import scipy.signal as sig
from determine_domain import get_ilat_ilon
from glob import glob

ilat, ilon = get_ilat_ilon()
ilonstr = [str(i) for i in ilon]
ilatstr = [str(i) for i in ilat]
npoints = len(ilat)

fileout = '/home1/scratch/nbarrier/clim_oope.nc'
dsclim = xr.open_dataset(fileout)
dsclim = dsclim['clim'].values
dsclim = dsclim[:, ilat, ilon, :]
print(dsclim.shape)
    
dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/debugged_corr_mask/output/'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

filelist = np.sort(glob('/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/debugged_corr_mask/output/*OOPE*nc'))
ntime = len(filelist) * 12

output = np.zeros((100, npoints, ntime), dtype=np.float)

index = np.arange(12)

for f in filelist:

    print('Processing file ', f)

    data = xr.open_dataset(f)
    data = data.isel(community=0)
    data = data['OOPE'].values[:, ilat, ilon, :]
    output[:, :, index] = (data - dsclim).T
    index += 12

fileout = '/home1/scratch/nbarrier/anoms_oope.nc'
dsout = xr.Dataset()
dsout['anom'] = (['size', 'points', 'time'], output)
#dsout.attrs['ilon'] = ', '.join(ilonstr)
#dsout.attrs['ilat'] = ', '.join(ilatstr)
dsout.to_netcdf(fileout)
