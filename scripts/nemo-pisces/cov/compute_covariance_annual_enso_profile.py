import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

index = xr.open_dataset('%s/yearly_nino_index_1950_2019.nc' %dirin)
enso = index['nino'].values
ensoy = index['year'].values

cdftime = utime("seconds since 1900-01-01 00:00:00", "noleap")

varlist = ['GOC',
'O2',
'PAR',
'PHY2',
'PHY',
'so',
'thetao',
'uocetr_eff',
'vocetr_eff',
'ZOO2',
'ZOO', 'CHL', 'PLK']

varlist = ['uo', 'vo']


for varname in varlist:

    print(varname)

    data = xr.open_dataset('%s/detrended_profile_yearly_%s.nc' %(dirin, varname))
    year = data['year'].values
    lon = data['lon'].values
    depth = data['depth'].values

    ymin = np.max([ensoy.min(), year.min()])
    ymax = np.min([ensoy.max(), year.max()])

    iok_enso = np.nonzero((ensoy <= ymax) & (ensoy >= ymin))[0]
    iok = np.nonzero((year <= ymax) & (year >= ymin))[0]
    
    data = data.isel(year=iok)
    enso = enso[iok_enso]
    ensoy = ensoy[iok_enso]
    N = len(enso)

    dimnames = data[varname].dims

    data = data[varname].to_masked_array()

    isea = np.nonzero(data.mask[0] == False)  
    cov = np.zeros(data.shape[1:])

    # data =  time, com, size, lat, lon
    for i, j in zip(isea[0], isea[1]):
        cov[i, j] = np.cov(data[:, i, j], enso, ddof=1)[0, 1]

    fileout = '%s/covariance_yearly_enso_profile_%s.nc' %(dirin, varname)
    output = xr.Dataset()
    output['covariance'] = (dimnames[1:], cov)
    output['lon'] = (['lon'], lon)
    output['depth'] = (['depth'], depth)
    output.attrs['file'] = os.path.realpath(__file__)
    output.attrs['date'] = str(datetime.datetime.today())
    output.to_netcdf(fileout)
