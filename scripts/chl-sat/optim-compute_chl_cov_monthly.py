import sys
sys.path.append("/home1/datahome/nbarrier/apecosm/configs/apecosm_orca1/diags/nino")
import os.path
import numpy as np
import xarray as xr
import logging
from extract_nino import read_index
from glob import glob

# Overlapping index for data (index) and nino index (nindex)
index = np.arange(0, 267 + 1)
nindex = np.arange(572, 839 + 1)

# read the nino index
daten, nino = read_index()
daten = daten[nindex]
nino = nino[nindex]

# reads the clim to extract data mask
dirin = '/home1/datawork/nbarrier/chl-data/'
dataclim = xr.open_dataset("%s/ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-CLIM-fv5.0.nc" %dirin)
nlat = dataclim.dims['lat']
nlon = dataclim.dims['lon']
lat = dataclim['lat'].values
lon = dataclim['lon'].values

filelist = np.sort(glob('%s/anom_*nc' %dirin))
filelist = filelist[index]

output = np.zeros((nlat, nlon), dtype=np.float)
counter = np.zeros((nlat, nlon), dtype=np.int)

N = len(filelist)
for i in range(N):
    f = filelist[i]
    print('Processing file %s' %f, daten[i])
    data = xr.open_dataset(f)
    data = np.squeeze(data['chlor_a'].values)
    data *= nino[i]
    
    iok = np.nonzero(np.logical_not(np.isnan(data)))

    output[iok] += data[iok]
    counter[iok] += 1

output = output / (counter - 1)
output = np.ma.masked_where(counter - 1 <= 0, output)

dsout = xr.Dataset()
dsout['covar'] = (['lat', 'lon'], output)
dsout['counter'] = (['lat', 'lon'], counter)
dsout['lat'] = (['lat',], lat)
dsout['lon'] = (['lon',], lon)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('%s/covariance_chl_nino_monthly_obs.nc' %dirin)
