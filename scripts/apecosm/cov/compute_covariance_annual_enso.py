import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime
    
dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/corr_mask/cyc4/output/yearly/apecosm/'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'
dirin = dirout

prefix = 'debugged_corr_mask'

index = xr.open_dataset('%s/yearly_nino_index_1950_2019.nc' %dirin)
ensof = index['nino'].values
ensoy = index['year'].values

cdftime = utime("seconds since 1900-01-01 00:00:00", "noleap")

for varname in ['OOPE']:

    print("@@@@@@@@@@@@@@@@@@ ", varname)

    data = xr.open_dataset('%s/%s_detrended_global_annual_%s.nc' %(dirin, prefix, varname))
    print(data)
    year = data['time'].values
    print(year)

    ymin = np.max([ensoy.min(), year.min()])
    ymax = np.min([ensoy.max(), year.max()])

    print(ymin, ymax)

    iok_enso = np.nonzero((ensoy <= ymax) & (ensoy >= ymin))[0]
    iok = np.nonzero((year <= ymax) & (year >= ymin))[0]

    data = data.isel(time=iok)
    enso = ensof[iok_enso]
    N = len(enso)
    
    dimnames = data[varname].dims

    data = data[varname].to_masked_array()
    print(data.shape)
    print(enso.shape)

    # extract the future size of cov. output
    shape1 = data.shape[1:]
    ntime = data.shape[0]
    shapetot = np.prod(shape1)

    # resize data into (ntime, nall) arrays
    data = np.ma.reshape(data, (ntime, shapetot))
    isea = np.nonzero(data.mask[0] == False)  

    # extract the covariance into 1D form
    cov = np.zeros(data[0].shape)

    # data =  time, com, size, lat, lon
    for i in isea[0]:
        cov[i] = np.cov(data[:, i], enso, ddof=1)[0, 1]

    # reshape covariance into it's proper form
    cov = np.reshape(cov, (shape1))

    fileout = '%s/%s_covariance_yearly_enso_%s.nc' %(dirin, prefix, varname)
    output = xr.Dataset()
    output['covariance'] = (dimnames[1:], cov)
    output.attrs['file'] = os.path.realpath(__file__)
    output.attrs['date'] = str(datetime.datetime.today())
    output.to_netcdf(fileout)
