import xarray as xr
import numpy as np
import apecosm.misc as misc
import matplotlib.pyplot as plt
import sys
sys.path.append('../nino')
from extract_nino import read_index
import os.path
from datetime import datetime

pref = 'corr_mask'
pref = 'debugged_corr_mask'
dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

daten, nino = read_index()
nino = np.ma.masked_where(np.isnan(nino), nino)

thres = 1

inino = np.nonzero(nino > 1.75)[0]
inina = np.nonzero(nino < -1)[0]
ineutral = np.nonzero(np.abs(nino) < 0.5)[0]

ynino = daten // 100
mnino = (daten - 100 * ynino).astype(np.int)

yfinal = []

for y in np.unique(ynino):

    test = (ynino == y)
    test = test & (mnino >= 10) & (mnino <= 12)
    iok = np.nonzero(test)[0]
    temp = np.mean(nino[iok])
    if(temp > 1):
        yfinal.append(y)
yfinal.remove(1986)
yfinal = np.array(yfinal)

yfinal = yfinal[yfinal >= 1958]

#for varname in ['mort_day', 'OOPE', 'starvation', 'repfonct_day']:
for varname in ['OOPE']:
    
    filename = '%s/%s_%s_meridional_mean_anoms.nc' %(dirin, pref, varname)
    print(filename)
    data = xr.open_dataset(filename)
    lon = data['x'].values
    time = data['time'].values
    year = data['time.year'].values
    month = data['time.month'].values
    oope = data[varname].values
    ntime, nx, ncom, nw = oope.shape
    date = year * 100 + month

    nyears = 12 * 4

    output = []

    for y in yfinal:

        istart = np.nonzero((year == y) & (month == 1))[0][0]
        iend = np.nonzero((year == y + 3) & (month == 12))[0][0]
        iend += 1
        print(date[istart:iend])
        output.append(oope[istart:iend])

    output = np.array(output)
    print(output.shape)
    compo_mean = np.mean(output, axis=0)
    compo_std = np.std(output, axis=0)

    dsout = xr.Dataset()
    dsout['compo_mean'] = (['time', 'lon', 'com', 'size'], compo_mean)
    dsout['compo_std'] = (['time', 'lon', 'com', 'size'], compo_std)
    dsout['compo_mean'].attrs['compo_years'] = yfinal
    dsout['lon'] = (['lon'], lon)
    dsout['time'] = (['time'], np.arange(nyears))
    dsout.attrs['file'] = os.path.realpath(__file__)
    dsout.attrs['date'] = str(datetime.today())
    dsout.to_netcdf('%s/%s_hovmoller_composite_%s.nc' %(dirin, pref, varname))

