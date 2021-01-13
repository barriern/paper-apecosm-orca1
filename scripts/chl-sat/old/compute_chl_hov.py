import os.path
from glob import glob
import numpy as np
import xarray as xr

output = np.zeros((12, 4320, 8640), dtype=np.float)
counter = np.zeros((12, 4320, 8640), dtype=np.int)

dirin = '/home1/datawork/nbarrier/chl-data/'
filelist = np.sort(glob('%s/ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-[0-9]*.nc' %(dirin)))
print(filelist)

for f in filelist:
    data = xr.open_dataset(f)
    lat = data['lat'].values
    ilat = np.nonzero(np.abs(lat) <= 2)[0]
    data = data.isel(lat=ilat)
    data = data.mean(dim='lat', skipna=True)
    data.attrs['file'] = os.path.realpath(__file__)

    fout = os.path.basename(f)
    fout = '%s/hov_%s' %(dirin, fout)

    data.to_netcdf(fout)
 
