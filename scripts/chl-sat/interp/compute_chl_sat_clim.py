import os.path
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from glob import glob

dirin = 'data/'
nlat = 180
nlon = 360
nmonths = 12

clim = np.zeros((nmonths, nlat, nlon), dtype=np.float)
counter = np.zeros((nmonths, nlat, nlon), dtype=np.int)

for y in range(1998, 2020):
    filelist = np.sort(glob('%s/*-%d*nc' %(dirin, y)))

    data = xr.open_mfdataset(filelist, combine='by_coords')
    lon = data['lon'].values
    lat = data['lat'].values
    count = data['count_a'].to_masked_array()
    data = data['chlor_a'].to_masked_array()
    
    data = np.ma.masked_where(count < (100 / 3), data)

    iok = np.nonzero(1 - np.ma.getmaskarray(data))  # iok = index where data is ocean
    clim[iok] += data[iok]
    counter[iok] += 1

clim = clim / counter
clim = np.ma.masked_where(counter==0, clim)

dsout = xr.Dataset()
dsout['clim_chl'] = (['month', 'lat', 'lon'], clim)
dsout['counter'] = (['month', 'lat', 'lon'], counter)
dsout['lon'] = (['lon',], lon)
dsout['lat'] = (['lat',], lat)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('%s/clim_chl_monthly_obs.nc' %dirin)
