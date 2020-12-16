import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime
from remove_trend_yearly import detrend

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'
dirin = 'data/'

index = xr.open_dataset('%s/yearly_nino_index_1950_2019.nc' %dirin)
index = index.isel(year=slice(0, -1))
enso = index['nino'].values
ensoy = index['year'].values

cdftime = utime("seconds since 1900-01-01 00:00:00", "noleap")

varlist = ['sst']

for varname in varlist:

    print(varname)

    data = xr.open_dataset('%s/yearly_mean_%s.nc' %(dirin, varname))
    year = data['time'].values
    lon = data['longitude'].values
    lat = data['latitude'].values

    ymin = 1958
    ymax = 2017

    iok_enso = np.nonzero((ensoy <= ymax) & (ensoy >= ymin))[0]
    iok = np.nonzero((year <= ymax) & (year >= ymin))[0]
    year = year[iok]
    
    data = data.isel(time=iok)

    enso = enso[iok_enso]
    enso = np.ma.masked_where(np.isnan(enso), enso)
    ensoy = ensoy[iok_enso]
    N = len(enso)

    dimnames = data[varname].dims

    data = data[varname].to_masked_array()

    print(data.shape, year.shape, year)
    print(enso.shape, ensoy.shape, ensoy)

    data = detrend(year, data) 

    isea = np.nonzero(data.mask[0] == False)  
    cov = np.zeros(data.shape[1:])

    # data =  time, com, size, lat, lon
    for i, j in zip(isea[0], isea[1]):
        cov[i, j] = np.cov(data[:, i, j], enso, ddof=1)[0, 1]

    fileout = '%s/covariance_yearly_enso_surface_%s.nc' %(dirin, varname)
    output = xr.Dataset()
    output['covariance'] = (dimnames[1:], cov)
    output.coords['longitude'] = lon
    output.coords['latitude'] = lat
    output.attrs['file'] = os.path.realpath(__file__)
    output.attrs['date'] = str(datetime.datetime.today())
    output.to_netcdf(fileout)
