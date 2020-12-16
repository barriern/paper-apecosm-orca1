import numpy as np
from datetime import date
import os.path
import numpy as np
from eofs.standard import Eof
import matplotlib.pyplot as plt
import cartopy.crs as crs
from scipy import stats
import xarray as xr
from remove_trend_yearly import detrend
import apecosm.ts as ts
import os.path
from datetime import datetime
import sys

sys.path.append('../nino')
from compute_composites import get_compo_year

yfinal = get_compo_year()
print(yfinal)


dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'
zmax = 200

varlist = ['thetao', 'ZOO', 'ZOO2', 'PHY2', 'O2', 'GOC', 'PAR', 'PHY', 'so']

for varname in varlist:

    print(varname)

    filein = '%s/hovmoller_data_0-%d_%s_anoms.nc' %(dirout, zmax, varname)
    data = xr.open_dataset(filein)

    lon = data['x'].values
    time = data['year'].values
    year = data['year.year'].values
    month = data['year.month'].values
    oope = data[varname].values
    ntime, nx = oope.shape
    date = year * 100 + month

    nyears = 12 * 4

    output = []

    for y in yfinal:

        istart = np.nonzero((year == y) & (month == 1))[0][0]
        iend = np.nonzero((year == y + 3) & (month == 12))[0][0]
        iend += 1
        output.append(oope[istart:iend])

    output = np.array(output)

    compo_mean = np.mean(output, axis=0)
    compo_std = np.std(output, axis=0)

    dsout = xr.Dataset()
    dsout['compo_mean'] = (['time', 'lon'], compo_mean)
    dsout['compo_std'] = (['time', 'lon'], compo_std)
    dsout['compo_mean'].attrs['compo_years'] = yfinal
    dsout['lon'] = (['lon'], lon)
    dsout['time'] = (['time'], np.arange(nyears))
    dsout.attrs['file'] = os.path.realpath(__file__)
    dsout.attrs['date'] = str(datetime.today())
    dsout.to_netcdf('%s/hovmoller_composite_0-%d_%s.nc' %(dirout, zmax, varname))
