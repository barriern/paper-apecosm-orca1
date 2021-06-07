import xarray as xr
from glob import glob
import numpy as np
import os.path

# extract month_duration (in seconds) in the format (12, 1, 1, 1)
month_duration = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
month_duration = month_duration[:, np.newaxis, np.newaxis, np.newaxis]  # time, lat, lon, w
month_duration *= 24 * 60 * 60

# open mesh file to extract x,y dimensions
mesh = xr.open_dataset('../../../data/mesh_mask_eORCA1_v2.2.nc')
ny = mesh.dims['y']
nx = mesh.dims['x']

dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data/integrated'

# list of variables to process
varnames = ['madv_trend']

for v in varnames:
    print(v)

    # extract the list of files
    pattern = '%s/*%s*nc' %(dirin, v)
    filelist = np.sort(glob(pattern))
    
    # extract values of previous last step (init. 0)
    previous = np.zeros((ny, nx, 100))
    
    clim = xr.open_dataset('%s/clim_%s.nc' %(dirout, v))
    clim = clim[v].to_masked_array()

    # filelist = filelist[:2]

    # loop over all the files
    for f in filelist:
        print(f)
        data = xr.open_dataset(f)
        time = data['time']
        # extract first community
        data = data.isel(community=0)
        data = data[v].to_masked_array()
        data = data - clim
       
        # computes cumulated sum over time
        # add values from the previous time step
        data = np.cumsum(data * month_duration, axis=0) + previous

        # Save the data in a given directory
        fout =  os.path.basename(f)
        fout = '%s/integrated_%s' %(dirout, fout)

        dsout = xr.Dataset()
        dsout['time'] = time
        dsout[v] = (['time', 'y', 'x', 'w'], data)
        dsout.attrs['file'] = os.path.realpath(__file__)
        dsout.to_netcdf(fout, unlimited_dims='time')
        
        # increment the previous value
        previous = data[-1]
