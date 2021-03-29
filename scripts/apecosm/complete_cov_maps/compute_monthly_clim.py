import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime
import scipy.signal as sig
from determine_domain import get_ilat_ilon
from glob import glob

ilat, ilon = get_ilat_ilon()

import sys
sys.path.append('/home1/datahome/nbarrier/apecosm/configs/apecosm_orca1/diags/nino/')
import extract_nino
    
dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/debugged_corr_mask/output/'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

filelist = np.sort(glob('/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/debugged_corr_mask/output/*OOPE*nc'))

for f in filelist:

    print('Processing file ', f)

    data = xr.open_dataset(f)
    data = data.isel(community=0)
    data = data['OOPE']

    if f == filelist[0]:
        clim = data.values
    else:
        clim += data.values

clim /= len(filelist)
print(clim.shape)

fileout = '/home1/scratch/nbarrier/clim_oope.nc'
dsout = xr.Dataset()
dsout['clim'] = (['month', 'lat', 'lon', 'size'], clim)
dsout.to_netcdf(fileout)
