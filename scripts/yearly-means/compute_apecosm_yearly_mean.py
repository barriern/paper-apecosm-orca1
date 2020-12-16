''' Computes the yearly averages of Apecosm outputs. 
    Yearly averages are computed from May (month 05) to April (month 04).
    
    One file per year is saved, since Apecosm outputs are too heavy.
'''

import xarray as xr
import os.path
from datetime import datetime
import numpy as np
import subprocess


''' Computes the yearly mean for Apecosm outputs. One file per year is saved. '''
def compute_yearly_mean(year, month, value, varname):

    # extracts the year to process, remove the last one 
    # since incomplete mean
    ytt = np.unique(year)
    ytt = ytt[:-1]
    
    date = year * 100 + month

    nyears = len(ytt)

    dimnames = value.dims

    # Loop over the years to process
    for y in ytt:
        print(y)

        # extracts the time index to use in the averaging
        test = (year == y) & (month >=5)
        test = test | ((year == (y + 1)) & (month < 5))
        iok = np.nonzero(test)[0]
        
        # computes the time mean, stored as DataArray
        # output is now of size (y, x, com, size)
        output = value.isel(time=iok).mean(dim='time').values
        
        # create the output dataset
        dataout = xr.Dataset()
        dataout[varname] = (dimnames, output[np.newaxis])  # add shadow time dimension
        dataout[dimnames[0]] = ([dimnames[0]], [y]) 
        dataout.attrs['file'] = os.path.realpath(__file__)
        dataout.attrs['date'] = str(datetime.today())
        dataout.attrs['description'] = 'mean over %s' %(str(date[iok]))
        fileout = '%s/%s_yearly_mean_%s_year_%.4d.nc' %(dirout, prefix, varname, y)
        dataout.to_netcdf(fileout)

    return dataout

# list of apecosm simulation to process
preflist = ['debugged_corr_mask']

# list of variables to process
varlist = ['OOPE']

# loop over all simulations
for prefix in preflist:
    
    # extract input dir.
    dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/%s/output' %prefix

    # creates the output directory (input_dir/yearly)
    dirout = '%s/yearly' %dirin
    subprocess.run(["mkdir", "-p", dirout])
    print(dirout)

    # loop over al variables
    for var in varlist:

        print('++++++++++++++++++++++++++++++++++ Computing %s yearly mean' %var)
        # open all files in one raw and compute average
        data = xr.open_mfdataset('%s/*%s*nc' %(dirin, var), combine='by_coords')
        data = data[var]
        year = data['time.year'].values
        month = data['time.month'].values
        output = compute_yearly_mean(year, month, data, var)
