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
import scipy.signal as sig

''' Removes the linear trend of Y, considering x as time and y of dimensions (time, ...) '''
def detrend(y):

    return sig.detrend(y.T).T


if __name__ == '__main__':

    preflist = ['final-runs']

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

        '''
        # extracts the 3cm, 20cm and 90cm size classes
        bins = np.array([3, 20, 90]).astype(np.float) * 1e-2
        ilength = []
        for l in bins:
            dist = (length - l)**2
            ilength.append(np.argmin(dist))
        ilength = np.array(ilength)
        '''

        # loop over the variables
        #for varname in ['u_active', 'u_passive', 'v_active', 'v_passive']:
        #for varname in ['madv_trend', 'gamma1', 'repfonct_day', 'mort_day', 'mdiff_trend', 'zdiff_trnd', 'zadv_trend']:
        for varname in ['zdiff_trend', 'zadv_trend', 'u_active', 'u_passive', 'v_active', 'v_passive']:
        #for varname in ['OOPE']:

            print('@@@@@@@@@@@@@ ', varname)
            pattern = '%s/*%s*_year_*.nc' %(dirin, varname)

            # opens the yearly OOPE file and extracts the size classes
            data = xr.open_mfdataset(pattern, combine='by_coords')
            data = data.isel(community=0) # time, y, x, length
            year = data['time'].values
            
            # extracts the dimension names
            dimnames = data[varname].dims
            
            # extracts the OOPE and detrend it using year
            dens = data[varname].to_masked_array()
            year = np.unique(year)

            ilat, ilon = np.nonzero(dens[0, :, :, 0].mask == False)
            dens[:, ilat, ilon, :] = detrend(dens[:, ilat, ilon, :])

            # Save detrended variables into file
            dataout = xr.Dataset()
            dataout['time'] = (['time'], year)
            dataout[varname] = (dimnames, dens)
            #dataout['length'] = (['length'], length)

            dataout.attrs['script'] = os.path.realpath(__file__)
            dataout.attrs['date'] = str(date.today())

            dataout.to_netcdf('%s/%s_detrended_global_annual_epis_%s.nc' %(dirout, pref, varname))
