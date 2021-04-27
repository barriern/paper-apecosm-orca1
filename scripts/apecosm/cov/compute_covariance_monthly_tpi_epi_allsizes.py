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

index = xr.open_dataset('../../data/filt_tpi.nc')
ensof = index['tpi_filt'].to_masked_array()
ensoy = index['time'].to_masked_array()

print(ensof.shape, ensoy.shape)
iok = np.nonzero(ensof.mask == False)[0]

ensoy = ensoy[iok]
ensof = ensof[iok]
print(ensoy)

cdftime = utime("seconds since 1900-01-01 00:00:00", "noleap")

for varname in ['OOPE']:

    print("@@@@@@@@@@@@@@@@@@ ", varname)

    data = xr.open_mfdataset('%s/*nc' %(dirin))
    year = data['time.year'].values * 100 + data['time.month'].values 

    ymin = np.max([ensoy.min(), year.min()])
    ymax = np.min([ensoy.max(), year.max()])

    iok_enso = np.nonzero((ensoy <= ymax) & (ensoy >= ymin))[0]
    iok = np.nonzero((year <= ymax) & (year >= ymin))[0]

    print(ymin, ymax)

    data = data.isel(time=iok)
    enso = ensof[iok_enso]
    N = len(enso)
    
    dimnames = data[varname].dims  # time, y, x, w
    dimnames = dimnames.T  # w, x, y, time

    data = data[varname].to_masked_array()  # time, y, x, w
    tmask = data.mask
    clim, data = ts.get_monthly_clim(data)
    del(clim)
    data = np.ma.masked_where(tmask, data)
    del(tmask)
    data = np.ma.transpose(data)  # w, x, y, time
    print(data.shape)
    print(enso.shape)

    # extract the future size of cov. output
    shape1 = data.shape[:-1]  # w, x, y
    ntime = data.shape[-1]  
    shapetot = np.prod(shape1)

    # resize data into (ntime, nall) arrays
    data = np.ma.reshape(data, (shapetot, ntime))
    isea = np.nonzero(data.mask[:, 0] == False)  

    # detrending of time-series
    data[isea, :] = sig.detrend(data[isea, :])

    # extract the covariance into 1D form
    cov = np.zeros(data[:, 0].shape)

    # data =  time, com, size, lat, lon
    for i in isea[0]:
        cov[i] = np.cov(data[i, :], enso, ddof=1)[0, 1]

    # reshape covariance into it's proper form
    cov = np.reshape(cov, (shape1))

    fileout = '%s/%s_covariance_monthly_tpi_epis_%s.nc' %(dirin, prefix, varname)
    output = xr.Dataset()
    output['covariance'] = (dimnames[:-1], cov)
    output.attrs['file'] = os.path.realpath(__file__)
    output.attrs['date'] = str(datetime.datetime.today())
    output.to_netcdf(fileout)
