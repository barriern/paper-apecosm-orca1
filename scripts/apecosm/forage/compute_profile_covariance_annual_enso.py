import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime
import scipy.signal as sig

# reading enso index
dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

index = xr.open_dataset('%s/yearly_nino_index_1950_2019.nc' %dirin)
ensof = index['nino'].values
ensoy = index['year'].values

cdftime = utime("seconds since 1900-01-01 00:00:00", "noleap")

data = xr.open_mfdataset('%s/yearly_mean_forage_year_*nc' %dirin, combine='by_coords')
year = data['time'].values

ymin = np.max([ensoy.min(), year.min()])
ymax = np.min([ensoy.max(), year.max()])

iok_enso = np.nonzero((ensoy <= ymax) & (ensoy >= ymin))[0]
iok = np.nonzero((year <= ymax) & (year >= ymin))[0]

data = data.isel(time=iok)
enso = ensof[iok_enso]
N = len(enso)

dimnames = data['FORAGE'].dims

forage = data['FORAGE'].to_masked_array()

forage = forage.T
ni, nj, nk, nz, nw, nt = forage.shape
i, j, k, z, w, t  = np.nonzero(forage.mask == False)

cov = np.zeros((ni, nj, nk, nz, nw))

for temp in zip(i, j, k, z, w):

    ii, jj, kk, zz, ww = temp
    
    ts = sig.detrend(forage[ii, jj, kk, zz, ww])
    cov[ii, jj, kk, zz, ww] = np.cov(ts, enso, ddof=1)[0, 1]

fileout = '%s/profile_covariance_yearly_enso.nc' %dirin
output = xr.Dataset()
output['covariance'] = (dimnames[1:], cov.T)
output.attrs['file'] = os.path.realpath(__file__)
output.attrs['date'] = str(datetime.datetime.today())
output.to_netcdf(fileout)
