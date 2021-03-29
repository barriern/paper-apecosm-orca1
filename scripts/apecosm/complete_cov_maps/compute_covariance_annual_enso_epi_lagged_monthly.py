import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime
import scipy.signal as sig
from determine_domain import get_ilat_ilon

import sys
sys.path.append('/home1/datahome/nbarrier/apecosm/configs/apecosm_orca1/diags/nino/')
import extract_nino

datenino, nino = extract_nino.read_index()
iok = np.nonzero((datenino <= 201812) & (datenino >= 195801))
#iok = np.nonzero((datenino <= 195912) & (datenino >= 195801))
datenino = datenino[iok]
nino = nino[iok]
    
dirout = '/home1/scratch/nbarrier/'

prefix = 'debugged_corr_mask'

cdftime = utime("seconds since 1900-01-01 00:00:00", "noleap")

for varname in ['OOPE']:

    print("@@@@@@@@@@@@@@@@@@ ", varname)

    filename = '/home1/scratch/nbarrier/anoms_oope.nc'
    data = xr.open_dataset(filename)
    data = data['anom']
    print(data.shape)

    nsize, npoints, ntime = data.shape
    
    nlags = 2 * (ntime - 1) + 1
    lags = np.arange(nlags) - (ntime - 1)
    ilags = np.nonzero(np.abs(lags) <= 5 * 12)[0]
    outlags = lags[ilags]
    nlagout = len(outlags)

    output = np.zeros((nsize, npoints, nlagout), dtype=np.int)

    for s in range(nsize):
        print('%s / %d' %(s, nsize))
        for j in range(npoints):
            y = data.isel(size=s, points=j)
            y = sig.detrend(y)
            output[s, j, :] = sig.correlate(y, nino)[ilags]

    print('End Loop')

    dsout = xr.Dataset()
    dsout['cov'] = (['w', 'points', 'lags'], output)
    dsout['lags'] = (['lags'], outlags)
    encoding = dict()
    encoding['cov'] = {} # {"zlib": True, "complevel": 9}
    dsout.to_netcdf('%s/lagged_covariance_monthly_epi_allsizes.nc' %(dirout), encoding=encoding)

    print('End program')
