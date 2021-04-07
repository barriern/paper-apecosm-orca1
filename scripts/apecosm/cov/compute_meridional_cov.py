import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime
import sys
import scipy.signal as sig
import re
from glob import glob
import apecosm.ts as ts

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'
lonmax = -150

pattern = '%s/zonal_mean_%.f_final-runs_(.*).nc' %(dirin, lonmax)
print(pattern)
regexp = re.compile(pattern)
filelist = glob('%s/*zonal_mean_%.f_final-runs*nc' %(dirin, lonmax))
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

for varname in varlist:
    
    print('Processing variable ', varname)

    filename = '%s/zonal_mean_%.f_final-runs_%s.nc' %(dirin, lonmax, varname)
    print(filename)
    data = xr.open_dataset(filename)
    data = data.isel(community=0)
    lat = data['y'].values
    oope = data[varname].to_masked_array()  # time, y, comm, w
    mask = np.ma.getmaskarray(oope)
    clim, oope = ts.get_monthly_clim(oope)
    oope = np.ma.masked_where(mask == True, oope)

    oope = oope.T  # w, lat, time
    nw, nlat, ntime = oope.shape

    iok = np.nonzero(np.ma.getmaskarray(oope[0, :, 0]) == False)[0]

    oope[:, iok, :] = sig.detrend(oope[:, iok, :])
    lags = sig.correlation_lags(ntime, nino.size)
    ilags = np.nonzero(np.abs(lags) <= 5*12)[0]
    lags = lags[ilags]
    nlags = len(lags)

    covariance = np.zeros((nw, nlat, nlags), dtype=np.float)

    for s in range(3):
        for i  in iok:
            temp = oope[s, i]
            covariance[s, i, :] = sig.correlate(temp, nino)[ilags]

    dsout = xr.Dataset()
    dsout['covariance'] = (['w', 'lat', 'lags'], covariance)
    dsout['lags'] = (['lags'], lags)
    dsout['lat'] = (['lat'], lat)
    dsout.to_netcdf('%s/meridional_covariance_%.f_%s.nc' %(dirin, lonmax, varname))
