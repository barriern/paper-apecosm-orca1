import pandas as pd
import numpy as np
import pylab as plt
import xarray as xr
import datetime
import os.path

def read_index(filename='index/oni.data', skipfooter=8, na=-99.90):

    data = pd.read_csv(filename, skiprows=1, skipfooter=skipfooter, engine='python', header=None, index_col=0, delim_whitespace=True, na_values=na)

    years = data.index
    months = np.arange(1, 13)

    mm, yy = np.meshgrid(months, years)
    date = np.ravel(yy*100 + mm)
    ts = np.ravel(data.values)

    year = date // 100

    return date, ts

if __name__ == '__main__':

    date, ts = read_index()
    year = (date / 100).astype(np.int)
    month = (date - year * 100).astype(np.int)

    output = []

    for y in np.unique(year):

        test = (year == (y)) & ((month == 11) | (month == 12))
        test = test | ((year == y + 1) & (month == 1))
        iok = np.nonzero(test)
        output.append(ts[iok].mean())


    ds = xr.Dataset()
    ds['nino'] = (['year'], np.array(output))
    ds['year'] = (['year'], np.unique(year))
    ds['nino'].attrs['description'] = 'Value(1950) = Average(1950-11, 1950-12, 1951-01)'
    ds.attrs['date'] = str(datetime.datetime.today())
    ds.attrs['file'] = os.path.realpath(__file__)
    ds.to_netcdf('./yearly_nino_index_%d_%d.nc' %(year.min(), year.max()))
