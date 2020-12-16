''' Computes the yearly averages of NEMO-Pisces outputs. 
    Yearly averages are computed from May (month 05) to April (month 04)
'''

import xarray as xr
import os.path
from datetime import datetime
import numpy as np

''' Computes the yearly mean. Returns a Dataset '''
def compute_yearly_mean(year, month, value, varname):

    # extracts the years to process
    ytt = np.unique(year)
    ytt = ytt[:-1]   # remove the last year (since not full average)
    nyears = len(ytt)
    
    # converts years/month into date
    date = year * 100 + month

    # get the size of the output array
    outdims = value.isel(time_counter=slice(0, nyears)).shape
    
    # extracts the names of the array dimension
    dimnames = value.dims

    # init the output array
    output = np.zeros(outdims)

    cpt = 0
    for y in ytt:
        
        # extracts the time index that corresponds to one year (from J to A)
        test = (year == y) & (month >=5)
        test = test | ((year == (y + 1)) & (month < 5))
        iok = np.nonzero(test)[0]
        # computes the mean over the time dimension
        output[cpt] = value.isel(time_counter=iok).mean(dim='time_counter')
        cpt += 1

    # converts numpy to dataset
    dataout = xr.Dataset()
    dataout[varname] = (dimnames, output)
    dataout[dimnames[0]] = ([dimnames[0]], ytt)

    return dataout


dirin = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/JRA_CO2/'

pref = 'ptrc_T'
varnames = ['GOC', 'O2', 'PHY', 'PHY2', 'POC', 'ZOO', 'ZOO2']

pref = 'grid_T'
varnames = ['thetao', 'so']

pref = 'grid_U'
varnames = ['uocetr_eff']

pref = 'grid_V'
varnames = ['vocetr_eff']

pref = 'diad_T'
varnames = ['PAR']

pref = 'add_T'
varnames = ['NO3', 'PO4', 'NH4', 'DCHL', 'NCHL']

pref = 'speed_U'
varnames = ['uo']

#pref = 'speed_V'
#varnames = ['vo']
    
# Reads the entire file that corresponds to the prefix (grid_T, ptrc_T, etc.)
dataglob = xr.open_mfdataset('%s/*%s*nc' %(dirin, pref), combine='by_coords')
year = dataglob['time_counter.year'].values
month = dataglob['time_counter.month'].values

# sets the output directory
dirout = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/corr_mask/output/yearly/pisces'

# loop over all variables
for v in varnames:

    # if file exists, nothing is done.
    fileout = '%s/yearly_mean_%s.nc' %(dirout, v)
    if(os.path.isfile(fileout)):
        continue

    print('++++++++++++++++++++++++++++++++++ Computing %s yearly mean' %v)
    
    # extracts the data as xr.DataArray
    data = dataglob[v]
    
    # computes yearly mean and return a xr.Dataset
    output = compute_yearly_mean(year, month, data, v)
    
    # add attributes + save in netcdf
    output.attrs['file'] = os.path.realpath(__file__)
    output.attrs['date'] = str(datetime.today())
    output.to_netcdf(fileout)
