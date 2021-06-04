import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime
import sys
import scipy.signal as sig
import re
from glob import glob

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

pattern = '%s/final-runs_(.*)_meridional_mean_anoms.nc' %dirin
regexp = re.compile(pattern)
filelist = glob('%s/*meridional_mean_anoms*nc' %dirin)
filelist.sort()
varlist = []
for f in filelist:
    varname = regexp.match(f).groups()[0]
    varlist.append(varname)

sys.path.append('/home1/datahome/nbarrier/apecosm/configs/apecosm_orca1/diags/paper-apecosm-orca1/scripts/nino')
from extract_nino import read_index

dmin, dmax = 195801, 201812
daten, nino = read_index(filename='/home1/datahome/nbarrier/apecosm/configs/apecosm_orca1/diags/nino/index/oni.data')
iok = np.nonzero((daten >= 195801) & (daten <= 201812))
daten, nino = daten[iok], nino[iok]
nnino = len(nino)

duration =  np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])  # days
duration = np.tile(duration, (61 , 1))
duration = np.ravel(duration)  # 732
duration = duration * 24 * 60 * 60 # h, m, s 

for varname in varlist:
    
    print('Processing variable ', varname)

    data = xr.open_dataset('%s/final-runs_%s_meridional_mean_anoms.nc' %(dirin, varname))
    data = data.isel(community=0)
    lon = data['x'].values
    oope = data[varname].to_masked_array()  # time, lon, w

    oope = np.cumcum(oope * duration, axis=0)

    oope = oope.T  # w, lon, time
    nw, nlon, ntime = oope.shape

    iok = np.nonzero(oope[0, :, 0].mask == False)[0]   # extract the lon where no NaNs

    oope[:, iok, :] = sig.detrend(oope[:, iok, :])  
    lags = sig.correlation_lags(ntime, nino.size)
    ilags = np.nonzero(np.abs(lags) <= 5*12)[0]
    lags = lags[ilags]
    nlags = len(lags)

    covariance = np.zeros((nw, nlon, nlags), dtype=np.float)

    for s in range(nw):
        for i  in iok:
            temp = oope[s, i]
            covariance[s, i, :] = sig.correlate(temp, nino)[ilags]

    dsout = xr.Dataset()
    dsout['covariance'] = (['w', 'lon', 'lags'], covariance)
    dsout['lags'] = (['lags'], lags)
    dsout['lon'] = (['lon'], lon)
    dsout.attrs['file'] = os.path.realpath(__file__)
    dsout.to_netcdf('%s/equatorial_covariance_%s.nc' %(dirin, varname))
