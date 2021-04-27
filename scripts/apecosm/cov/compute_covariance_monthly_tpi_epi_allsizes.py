import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime
import apecosm.ts as ts
import scipy.signal as sig
    
dirin = '/home1/scratch/nbarrier'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

prefix = 'final-runs'

# read the data of filtered TPI
index = xr.open_dataset('../../data/filt_tpi.nc')
ensof = index['tpi_filt'].to_masked_array()
ensoy = index['time'].to_masked_array()

# extracts the time series of non-masked data.
print(ensof.shape, ensoy.shape)
iok = np.nonzero(ensof.mask == False)[0]

ensoy = ensoy[iok]
ensof = ensof[iok]
print(ensoy)

cdftime = utime("seconds since 1900-01-01 00:00:00", "noleap")

for varname in ['OOPE']:

    print("@@@@@@@@@@@@@@@@@@ ", varname)

    data = xr.open_mfdataset('%s/temp*nc' %(dirin))
    year = data['time.year'].values * 100 + data['time.month'].values 

    ymin = np.max([ensoy.min(), year.min()])
    ymax = np.min([ensoy.max(), year.max()])

    # extract the time-series ovelap
    iok_enso = np.nonzero((ensoy <= ymax) & (ensoy >= ymin))[0]
    iok = np.nonzero((year <= ymax) & (year >= ymin))[0]

    print(ymin, ymax)

    # extraction of good time-series
    data = data.isel(time=iok)
    enso = ensof[iok_enso]
    N = len(enso)
    
    dimnames = list(data.dims)  # time, y, x, w
    print(dimnames)
    dimnames = dimnames[::-1]  # w, x, y, time

    data = data[varname].to_masked_array()  # time, y, x, w
    clim, data = ts.get_monthly_clim(data)
    del(clim)
    print(data.shape)
    data = np.transpose(data)
    print(data.shape)

    # extract the future size of cov. output
    shape1 = data.shape[:-1]  # w, x, y
    ntime = data.shape[-1]  
    shapetot = np.prod(shape1)

    nw, nx, ny = shape1

    # extract the covariance into 1D form
    cov = np.zeros(shape1)  # w, y, x
    print(cov.shape)
    print(data.shape)

    # data =  time, com, size, lat, lon
    for w in range(nw):
        for i in range(nx):
            for j in range(ny):
                if np.any(np.isnan(data[w, i, j :])):
                    continue
                cov[w, i, j] = np.cov(sig.detrend(data[w, i, j :]), enso, ddof=1)[0, 1]

    fileout = '%s/%s_covariance_monthly_tpi_epis_%s.nc' %(dirin, prefix, varname)
    output = xr.Dataset()
    output['covariance'] = (dimnames, cov)  # w, x, y
    output.attrs['file'] = os.path.realpath(__file__)
    output.attrs['date'] = str(datetime.datetime.today())
    output.to_netcdf(fileout)

