import os.path
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from glob import glob

dirin = '/home1/datawork/nbarrier/chl-data/'

nlon = 8640
nmonths = 12

clim = np.zeros((nmonths, nlon), dtype=np.float)
counter = np.zeros((nmonths, nlon), dtype=np.int)

for y in range(1998, 2020):
    filelist = np.sort(glob('%s/hov*-%d*nc' %(dirin, y)))
    print(y, len(filelist))
    print(filelist)

    data = xr.open_mfdataset(filelist, combine='by_coords')
    lon = data['lon'].values
    data = data['chlor_a'].to_masked_array()
    print(data.shape)
    print(clim.shape)

    iok = np.nonzero(1 - np.ma.getmaskarray(data))  # iok = index where data is ocean
    clim[iok] += data[iok]
    counter[iok] += 1

clim = clim / counter
clim = np.ma.masked_where(counter==0, clim)

dsout = xr.Dataset()
dsout['clim_chl'] = (['month', 'lon'], clim)
dsout['counter'] = (['month', 'lon'], counter)
dsout['lon'] = (['lon',], lon)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('%s/hov_clim_chl_monthly_obs.nc' %dirin)
