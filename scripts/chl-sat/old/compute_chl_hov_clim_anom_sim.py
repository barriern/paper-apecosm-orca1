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

dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

data = xr.open_dataset('%s/simulated_surface_chl_hov.nc' %dirout)
year = data['time_counter.year'].values
month = data['time_counter.month'].values
date = year * 100 + month
iok = np.nonzero(date >= 199709)[0]
data = data.isel(time_counter=iok)

lon = data['x']
lon = (lon + 360) % 360
data['x'] = lon

data = data.sortby(data['x'])
data = data.sel(x=slice(130, 300))
lon = data['x']

year = np.arange(1998, 2019)

yeardata = data['time_counter.year'].values
monthdata = data['time_counter.month'].values - 1

clim = np.zeros((12, len(lon)))
print(clim.shape)

cpt = 0
for y in year:

    iok = np.nonzero(yeardata == y)[0]
    clim += data.isel(time_counter=iok)['chl'].to_masked_array()
    cpt += 1

clim /= cpt

anom = data['chl'] - clim[monthdata]

anom.to_netcdf('%s/anom_chl_simulated.nc' %dirout)

ccc = 0.2

plt.figure()
cs = anom.plot()
cs.set_clim(-ccc, ccc)
plt.savefig('anom_chl_simulated.png')


