from glob import glob
import numpy as np
import xarray as xr

output = np.zeros((12, 4320, 8640), dtype=np.float)
counter = np.zeros((12, 4320, 8640), dtype=np.int)

dirin = '/home1/datawork/nbarrier/chl-data/'

for y in range(1998, 2020):

    print("++++++++++++++++++++++++ processing ", y)

    filelist = np.sort(glob('%s/*-%.4d*nc' %(dirin, y)))
    print(filelist)
    data = xr.open_mfdataset(filelist, combine='by_coords')
    lat = data['lat'].values
    lon = data['lon'].values
    data = data['chlor_a'].to_masked_array()

    output[data.mask == False] += data[data.mask == False]
    counter += (1 - data.mask.astype(np.int)) # 0 if masked

print("writting output")

output = output / counter
output = np.ma.masked_where(counter == 0, output)

dataout = (['month', 'lat', 'lon'], output)
datacount = (['month', 'lat', 'lon'], counter)
datalat = (['lat', ], lat)
datalon = (['lon', ], lon)
dataout = xr.Dataset({'clim_chlor_a' : dataout, 'lat':datalat, 'lon': datalon, 'counter': datacount})
dataout.to_netcdf("ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-CLIM-fv5.0.nc")
