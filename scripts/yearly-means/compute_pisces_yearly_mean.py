import xarray as xr
import os.path
from datetime import datetime
import numpy as np

def compute_yearly_mean(year, month, value, varname):

    ytt = np.unique(year)
    ytt = ytt[:-1]   # remove the last year (since not full average)
    nyears = len(ytt)
    
    date = year * 100 + month

    outdims = value.isel(time_counter=slice(0, nyears)).shape
    dimnames = value.dims

    output = np.zeros(outdims)

    cpt = 0
    for y in ytt:

        test = (year == y) & (month >=5)
        test = test | ((year == (y + 1)) & (month < 5))
        iok = np.nonzero(test)[0]
        print(date[iok])
        output[cpt] = value.isel(time_counter=iok).mean(dim='time_counter')
        cpt += 1

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
    
dataglob = xr.open_mfdataset('%s/*%s*nc' %(dirin, pref), combine='by_coords')
year = dataglob['time_counter.year'].values
month = dataglob['time_counter.month'].values

dirout = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/corr_mask/output/yearly/pisces'

for v in varnames:

    fileout = '%s/yearly_mean_%s.nc' %(dirout, v)
    if(os.path.isfile(fileout)):
        continue

    print('++++++++++++++++++++++++++++++++++ Computing %s yearly mean' %v)

    data = dataglob[v]
    output = compute_yearly_mean(year, month, data, v)
    output.attrs['file'] = os.path.realpath(__file__)
    output.attrs['date'] = str(datetime.today())
    output.to_netcdf(fileout)
