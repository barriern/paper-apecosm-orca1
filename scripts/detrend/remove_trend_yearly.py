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

def detrend(x, y):

    # compute the anomaly
    y = y - np.mean(y, axis=0, keepdims=1)

    # compute the shape of the rightmost dimensions
    yshape = y.shape
    ntime = y.shape[0]
    otherdims = y.shape[1:]

    alldims = np.prod(otherdims)
    
    # reshape onto (time, alldims)
    y1d = np.reshape(y, (ntime, alldims))

    # compute the trends by looping on the other dimensions
    trend = [stats.linregress(x, ytemp) for ytemp in y1d.T]

    # extracting the slopes and intercept
    # dimensions = ncells
    slopes = np.array([t[0] for t in trend]) 
    intercept = np.array([t[1] for t in trend])

    # x dimesions = ntime
    # computing linear trend
    output = slopes[np.newaxis, :] * x[:, np.newaxis] + intercept[np.newaxis, :]
    
    # move array to the right shape
    output = np.reshape(np.array(output), yshape)
    output = y - output

    return output

if __name__ == '__main__':


    preflist = ['corr_mask', 'climTemp']
    preflist = ['climPlk']

    for pref in preflist:
        
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ', pref)

        dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/%s/output/yearly/' %pref
        dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

        pattern = '%s/../*ConstantFields.nc' %dirin
        print(pattern)

        filemesh = glob('%s/../*ConstantFields.nc' %dirin)[0]

        const = xr.open_dataset(filemesh)
        length = const['length'].values

        bins = np.array([3, 20, 90]).astype(np.float) * 1e-2
        ilength = []
        for l in bins:
            dist = (length - l)**2
            ilength.append(np.argmin(dist))
        ilength = np.array(ilength)

        for varname in ['OOPE']:

            print('@@@@@@@@@@@@@ ', varname)

            #data = xr.open_mfdataset('data/density.nc', combine='by_coords')
            pattern = '%s/*%s*.nc' %(dirin, varname)
            print(np.sort(glob(pattern)))
            data = xr.open_mfdataset(pattern, combine='by_coords')
            data = data.isel(w=ilength)
            year = data['time'].values
            dimnames = data[varname].dims
            dens = data[varname].to_masked_array()
            year = np.unique(year)
            dens = detrend(year, dens)

            dataout = xr.Dataset()
            dataout['time'] = (['time'], year)
            dataout[varname] = (dimnames, dens)
            dataout['bins'] = (['bins'], bins)

            dataout.attrs['script'] = os.path.realpath(__file__)
            dataout.attrs['date'] = str(date.today())

            dataout.to_netcdf('%s/%s_detrended_global_annual_%s.nc' %(dirout, pref, varname))
