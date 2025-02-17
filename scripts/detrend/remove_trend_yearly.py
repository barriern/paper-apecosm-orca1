''' Scripts that removes linear trend from Apecosm yearly outputs '''

import numpy as np
from glob import glob
import xarray as xr
from datetime import date
import os.path
import numpy as np
from eofs.standard import Eof
import matplotlib.pyplot as plt
import cartopy.crs as crs
from scipy import stats

''' Removes the linear trend of Y, considering x as time and y of dimensions (time, ...) '''
def detrend(y):

    output = y.T   # w, c, x, y, t 
    print(output.shape)

    iok = np.nonzero(output.mask == False)
    output[iok].shape
    
    return output.T

if __name__ == '__main__':

    preflist = ['debugged_corr_mask']

    # loop over all the simulations
    for pref in preflist:
        
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ', pref)

        # reconstruct input dir, and output dor
        dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/%s/output/yearly/' %pref
        dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

        # reads the constant fields and extracts the length
        pattern = '%s/../*ConstantFields.nc' %dirin
        filemesh = glob('%s/../*ConstantFields.nc' %dirin)[0]
        const = xr.open_dataset(filemesh)
        length = const['length'].values

        # extracts the 3cm, 20cm and 90cm size classes
        bins = np.array([3, 20, 90]).astype(np.float) * 1e-2
        ilength = []
        for l in bins:
            dist = (length - l)**2
            ilength.append(np.argmin(dist))
        ilength = np.array(ilength)

        # loop over the variables
        for varname in ['OOPE']:

            print('@@@@@@@@@@@@@ ', varname)
            pattern = '%s/*%s*.nc' %(dirin, varname)
            # opens the yearly OOPE file and extracts the size classes
            data = xr.open_mfdataset(pattern, combine='by_coords')
            data = data.isel(w=ilength)  # time, y, x, com, w
            year = data['time'].values
            
            # extracts the dimension names
            dimnames = data[varname].dims
            
            # extracts the OOPE and detrend it using year
            dens = data[varname].to_masked_array()
            year = np.unique(year)
            dens = detrend(dens)

            # Save detrended variables into file
            dataout = xr.Dataset()
            dataout['time'] = (['time'], year)
            dataout[varname] = (dimnames, dens)
            dataout['bins'] = (['bins'], bins)

            dataout.attrs['script'] = os.path.realpath(__file__)
            dataout.attrs['date'] = str(date.today())

            dataout.to_netcdf('%s/%s_detrended_global_annual_%s.nc' %(dirout, pref, varname))
