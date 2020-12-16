import xarray as xr
import os.path
from datetime import datetime
import numpy as np

def compute_yearly_mean(year, month, value, varname):

    ytt = np.unique(year)
    nyears = len(ytt)
    
    date = year * 100 + month

    outdims = value.isel(time=slice(0, nyears)).shape
    dimnames = value.dims

    output = np.zeros(outdims)

    cpt = 0
    for y in ytt:

        test = (year == y) & (month >=5)
        test = test | ((year == (y + 1)) & (month < 5))
        iok = np.nonzero(test)[0]
        print(date[iok])
        if(len(iok) != 12):
            print(":@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        output[cpt] = value.isel(time=iok).mean(dim='time')
        cpt += 1

    dataout = xr.Dataset()
    dataout[varname] = (dimnames, output)
    dataout[dimnames[0]] = ([dimnames[0]], ytt)

    return dataout


dirin = 'data/'
    
dataglob = xr.open_mfdataset('%s/HadISST_sst.nc' %(dirin), combine='by_coords')
year = dataglob['time.year'].values
month = dataglob['time.month'].values
varnames = ['sst']
longitude = dataglob['longitude'].values
latitude = dataglob['latitude'].values
print("@@@@@ running")

for v in varnames:

    fileout = 'yearly_mean_%s.nc' %v
    if(os.path.isfile(fileout)):
        continue

    print('++++++++++++++++++++++++++++++++++ Computing %s yearly mean' %v)

    data = dataglob[v]
    output = compute_yearly_mean(year, month, data, v)
    output.attrs['file'] = os.path.realpath(__file__)
    output.attrs['date'] = str(datetime.today())
    output.coords['longitude'] = longitude
    output.coords['latitude'] = latitude
    output.to_netcdf(fileout)
