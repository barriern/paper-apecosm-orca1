import xarray as xr
from glob import glob
import numpy as np
import os.path

# open mesh file to extract x,y dimensions
mesh = xr.open_dataset('../../../data/mesh_mask_eORCA1_v2.2.nc')
lat = mesh['gphit'].values[0]
ilat, ilon = np.nonzero(np.abs(lat) <= 40)
ilat = slice(ilat.min(), ilat.max() + 1)
mesh = mesh.isel(y=ilat)
ny = mesh.dims['y']
nx = mesh.dims['x']

dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data/integrated'

# list of variables to process
varnames = ['madv_trend', 'zadv_trend', 'zdiff_trend']
varnames = ['mdiff_trend']

for v in varnames:
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ", v)
    
    clim = xr.open_dataset('%s/clim_%s.nc' %(dirout, v))
    clim = clim.isel(y=ilat)
    clim = clim[v].to_masked_array()

    # extract the list of files
    pattern = '%s/*%s*nc' %(dirin, v)
    filelist = np.sort(glob(pattern))
    
    # loop over all the files
    for f in filelist:
        print(f)
        data = xr.open_dataset(f)
        time = data['time']
        # extract first community
        data = data.isel(community=0, y=ilat)   # time, lat, lon, w
        data = data[v].to_masked_array() - clim
        
        fout =  os.path.basename(f)
        fout = '%s/anoms_%s' %(dirout, fout)

        dsout = xr.Dataset()
        dsout['time'] = time
        dsout[v] = (['time', 'y', 'x', 'w'], data)
        dsout.attrs['file'] = os.path.realpath(__file__)
        dsout.to_netcdf(fout, unlimited_dims='time')
