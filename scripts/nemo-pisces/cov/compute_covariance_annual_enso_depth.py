import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

index = xr.open_dataset('%s/yearly_nino_index_1950_2019.nc' %dirin)
enso = index['nino'].values
ensoy = index['year'].values
print(ensoy)

cdftime = utime("seconds since 1900-01-01 00:00:00", "noleap")

varlist = ['thetao']
varlist = ['thetao', 'ZOO', 'ZOO2', 'PHY2', 'O2', 'PHY', 'PAR', 'GOC']
varlist = ['uocetr_eff', 'vocetr_eff']
varlist = ['CHL']
varlist = ['uo', 'vo']
zmax = 200

for varname in varlist:

    data = xr.open_dataset('%s/detrended_0-%d_yearly_%s.nc' %(dirin, zmax, varname))
    year = data['year'].values

    ymin = np.max([ensoy.min(), year.min()])
    ymax = np.min([ensoy.max(), year.max()])

    iok_enso = np.nonzero((ensoy <= ymax) & (ensoy >= ymin))[0]
    iok = np.nonzero((year <= ymax) & (year >= ymin))[0]

    data = data.isel(year=iok)
    enso = enso[iok_enso]
    ensoy = ensoy[iok_enso]
    print(ensoy)
    N = len(enso)

    dimnames = data[varname].dims

    data = data[varname].to_masked_array()

    isea = np.nonzero(data.mask[0] == False)  
    cov = np.zeros(data.shape[1:])
    print(len(isea[0]))

    # data =  time, com, size, lat, lon
    for i, j in zip(isea[0], isea[1]):
        cov[i, j] = np.cov(data[:, i, j], enso, ddof=1)[0, 1]

    fileout = '%s/covariance_yearly_enso_0-%d_%s.nc' %(dirin, zmax, varname)
    output = xr.Dataset()
    output['covariance'] = (dimnames[1:], cov)
    output.attrs['file'] = os.path.realpath(__file__)
    output.attrs['date'] = str(datetime.datetime.today())
    output.to_netcdf(fileout)
