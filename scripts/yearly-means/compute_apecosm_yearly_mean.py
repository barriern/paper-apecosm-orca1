import xarray as xr
import os.path
from datetime import datetime
import numpy as np
import subprocess

def compute_yearly_mean(year, month, value, varname):

    ytt = np.unique(year)
    date = year * 100 + month

    ytt = ytt[:-1]

    nyears = len(ytt)

    dimnames = value.dims

    for y in ytt:
        print(y)

        test = (year == y) & (month >=5)
        test = test | ((year == (y + 1)) & (month < 5))
        iok = np.nonzero(test)[0]
        output = value.isel(time=iok).mean(dim='time').values
        print(output.shape)
    
        dataout = xr.Dataset()
        dataout[varname] = (dimnames, output[np.newaxis])
        dataout[dimnames[0]] = ([dimnames[0]], [y])
        dataout.attrs['file'] = os.path.realpath(__file__)
        dataout.attrs['date'] = str(datetime.today())
        dataout.attrs['description'] = 'mean over %s' %(str(date[iok]))
        fileout = '%s/%s_yearly_mean_%s_year_%.4d.nc' %(dirout, prefix, varname, y)
        dataout.to_netcdf(fileout)

    return dataout

preflist = ['debugged_corr_mask']

for prefix in preflist:

    dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/%s/output' %prefix

    dirout = '%s/yearly' %dirin
    subprocess.run(["mkdir", "-p", dirout])
    print(dirout)

    preflist = ['OOPE']

    for pref in preflist:

        print('++++++++++++++++++++++++++++++++++ Computing %s yearly mean' %pref)

        data = xr.open_mfdataset('%s/*%s*nc' %(dirin, pref), combine='by_coords')
        data = data[pref]
        year = data['time.year'].values
        month = data['time.month'].values
        output = compute_yearly_mean(year, month, data, pref)
