import xarray as xr
from glob import glob
import numpy as np
import os.path
import scipy.signal as sig

# extract month_duration (in seconds) in the format (12, 1, 1, 1)
month_duration = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
month_duration = month_duration[:, np.newaxis, np.newaxis, np.newaxis]  # time, lat, lon, w
month_duration *= 24 * 60 * 60

month_duration = np.tile(month_duration, (61, 1))
month_duration = np.ravel(month_duration)
month_duration = month_duration[:, np.newaxis, np.newaxis, np.newaxis]

# open mesh file to extract x,y dimensions
mesh = xr.open_dataset('../../../data/mesh_mask_eORCA1_v2.2.nc')
ny = mesh.dims['y']
nx = mesh.dims['x']

dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data/integrated'
dirin = dirout

# list of variables to process
varnames = ['madv_trend', 'zdiff_trend', 'mdiff_trend']
varnames = ['zdiff_trend']

for v in varnames:
    print(v)

    # extract the list of files
    pattern = '%s/anoms*%s*nc' %(dirin, v)
    print(pattern)
    filelist = np.sort(glob(pattern))

    data = xr.open_mfdataset(filelist)
    time = data['time']
    ntime = data.dims['time']
    dimnames = data[v].dims
    data = np.cumsum(data[v].to_masked_array() * month_duration, axis=0)   # time, lat, lon, w

    '''
    ilon, ilat = np.nonzero(np.ma.getmaskarray(data[0, :, :, 0]) == False)

    for w in range(100):
        for i, j in zip(ilon, ilat):
            data[w, i, j, :] = sig.detrend(data[w, i, j, :])
            data[w, i, j, :] = np.cumsum(data[w, i, j, :] * month_duration)
    '''

    # Save the data in a given directory
    #fout =  os.path.basename(f)
    fout = '%s/integrated_%s.nc' %(dirout, v)

    dsout = xr.Dataset()
    dsout['time'] = time
    dsout[v] = (dimnames, data)
    dsout.attrs['file'] = os.path.realpath(__file__)
    dsout.to_netcdf(fout, unlimited_dims='time')

    del(data)
