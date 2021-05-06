import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime
import scipy.signal as sig
import sys
sys.path.append('../../nino')
from extract_nino import read_index
from glob import glob

date, nino = read_index('../../data/index/oni.data')

# reading enso index
dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

filelist = glob('%s/detrended_equatorial_*_anomalies.nc' %dirin)
filelist.remove('%s/detrended_equatorial_forage_anomalies.nc' %dirin)

dataglob = xr.open_mfdataset(filelist, combine='by_coords')
dataglob['PLK_anoms'] = dataglob['GOC_anoms'] + dataglob['PHY2_anoms'] + dataglob['ZOO_anoms'] + dataglob['ZOO2_anoms']
dataglob['CHL_anoms'] = dataglob['NCHL_anoms'] + dataglob['DCHL_anoms']
print(dataglob)

for varname in ['PLK', 'CHL', 'O2']:

    print("Process variable ", varname)

    data = dataglob['%s_anoms' %varname]
    dimnames = data.dims  # x, z, time
    data = data.values

    ns, nc, ntime = data.shape
    cov = np.zeros((ns, nc), dtype=float)

    for i in range(ns):
        for j in range(nc):
            cov[i, j] = np.cov(data[i, j], nino, ddof=1)[0, 1]

    fileout = '%s/zonal_monthly_covariance_monthly_%s_enso.nc' %(dirin, varname)
    output = xr.Dataset()
    output['covariance'] = (dimnames[:-1], cov)
    output.attrs['file'] = os.path.realpath(__file__)
    output.attrs['date'] = str(datetime.datetime.today())
    output.to_netcdf(fileout)
