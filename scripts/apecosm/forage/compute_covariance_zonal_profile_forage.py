import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime
import scipy.signal as sig
import sys
sys.path.append('../../nino')
from extract_nino import read_index

date, nino = read_index('../../data/index/oni.data')
print(date.shape)
print(nino.shape)

# reading enso index
dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

data = xr.open_dataset('%s/detrended_equatorial_forage_anomalies.nc' %(dirin))
data = data['forage_anoms']
dimnames = data.dims
data = data.values

ns, nc, nz, nx, nd, ntime = data.shape
cov = np.zeros((ns, nc, nz, nx, nd), dtype=float)

for i in range(ns):
    for j in range(nc):
        for k in range(nz):
            for l in range(nx):
                for m in range(nd):
                    cov[i, j, k, l, m] = np.cov(data[i, j, k, l, m], nino, ddof=1)[0, 1]

fileout = '%s/zonal_monthly_covariance_yearly_enso.nc' %(dirin)
output = xr.Dataset()
output['covariance'] = (dimnames[:-1], cov)
output.attrs['file'] = os.path.realpath(__file__)
output.attrs['date'] = str(datetime.datetime.today())
output.to_netcdf(fileout)
