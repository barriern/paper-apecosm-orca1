import os.path
import xarray as xr
import numpy as np

data = xr.open_mfdataset("data/hovmoller_forage_year_*.nc", combine='by_coords')
year = data['time.year'].values
month = data['time.month'].values
date = year * 100 + month
forage = data['FORAGE']

ytt = np.unique(year)

for y in ytt:

    print(y)

    fout = 'data/forage_yearly_mean_%.4d.nc' %y

    test = (year == y) & (month >=5)
    test = test | ((year == (y + 1)) & (month < 5))
    iok = np.nonzero(test)[0]
    dtemp = date[iok]
    if(len(dtemp) != 12):
        print('Year %d discarded' %y)
        continue

    temp = forage.isel(time=iok).mean(dim='time', keepdims=True)
    temp['time'] = [y]
    temp.attrs['file'] = os.path.realpath(__file__)
    temp.to_netcdf(fout, unlimited_dims='time')
