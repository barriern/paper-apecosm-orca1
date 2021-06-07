import xarray as xr
from glob import glob
import numpy as np
import os.path

# open mesh file to extract x,y dimensions
mesh = xr.open_dataset('../../../data/mesh_mask_eORCA1_v2.2.nc')
ny = mesh.dims['y']
nx = mesh.dims['x']

dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data/integrated'

# list of variables to process
varnames = ['madv_trend', 'zadv_trend', 'zdiff_trend', 'zadv_trend']

for v in varnames:
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ", v)

    # extract the list of files
    pattern = '%s/*%s*nc' %(dirin, v)
    filelist = np.sort(glob(pattern))
    
    # extract values of previous last step (init. 0)
    output = np.zeros((12, ny, nx, 100))

    # loop over all the files
    for f in filelist:
        print(f)
        data = xr.open_dataset(f)
        time = data['time']
        # extract first community
        data = data.isel(community=0)   # time, lat, lon, comm, w
        data = data[v].to_masked_array()

        output += data

    output /= len(filelist)
       
    # Save the data in a given directory
    fout =  os.path.basename(f)
    fout = '%s/clim_%s.nc' %(dirout, v)

    dsout = xr.Dataset()
    dsout[v] = (['time', 'y', 'x', 'w'], output)
    dsout.attrs['file'] = os.path.realpath(__file__)
    dsout.to_netcdf(fout, unlimited_dims='time')
