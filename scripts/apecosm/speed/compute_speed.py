import xarray as xr
import numpy as np
from glob import glob
import os

dirin = "/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output"

ufiles = np.sort(glob('%s/*u_active*nc' %dirin))
vfiles = np.sort(glob('%s/*v_active*nc' %dirin))

nfiles = len(ufiles)
nfiles = 1

for p in range(nfiles):
    
    utemp = ufiles[p]
    vtemp = vfiles[p]

    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ", p)

    print("+++++ u file = ", os.path.basename(utemp))
    outputname = utemp.replace('u_active', 'mod_active')


    datau = xr.open_dataset(utemp)
    datau = datau['u_active']
    print('before ', datau.shape)
    datau = datau.isel(x=slice(1, None)) + datau.isel(x=slice(0, -1))
    datau = datau.isel(y=slice(1, None))
    print('after ', datau.shape)
    
    print("+++++ v file = ", os.path.basename(vtemp))
    datav = xr.open_dataset(vtemp)
    datav = datav['v_active']
    print('before ', datav.shape)
    datav = datav.isel(y=slice(1, None)) + datav.isel(y=slice(0, -1))
    datav = datav.isel(x=slice(1, None))
    print('after ', datav.shape)

    speed = (datau * datau + datav * datav)
    speed = speed.rename('speed')
    speed.to_netcdf(outputname, unlimited_dims='time')
    print(speed)


