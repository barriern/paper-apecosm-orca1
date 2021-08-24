import xarray as xr
from glob import glob
import numpy as np
import os.path

dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output/'
dirout = '/home1/scratch/nbarrier/'

latmax = 50

mesh = xr.open_dataset("../data/mesh_mask_eORCA1_v2.2.nc")
lon = mesh['glamt'].values[0]
lat = mesh['gphit'].values[0]

test = (np.abs(lat) <= latmax)
test2 = (lon >= 130) | (lon <= -60)
test = test & test2

ilat, ilon = np.nonzero(test == 1)
ilat = slice(ilat.min(), ilat.max() + 1)
ilon = slice(ilon.min(), ilon.max() + 1)

filelist = np.sort(glob('%s/*OOPE*nc' %dirin))
#filelist = filelist[:2]

print("ilat", ilat)
print("ilon", ilon)

for f in filelist:
    print("++++++++ processing ", f)

    bname = os.path.basename(f)

    data = xr.open_dataset(f)
    print(data['OOPE'].shape)
    #data = data.isel(community=3, x=ilon, y=ilat)['OOPE']
    data = data.isel(x=ilon, y=ilat)['OOPE']
    data.attrs['file'] = os.path.realpath(__file__)
    data.attrs['ilon'] = str(ilon)
    data.attrs['ilat'] = str(ilat)
    data.to_netcdf('%s/temp_%s' %(dirout, bname), unlimited_dims='time')
