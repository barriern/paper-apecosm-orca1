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
import misc as ts
import os.path
from datetime import datetime

dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'
zmax = 200

varlist = ['thetao', 'ZOO', 'ZOO2', 'PHY2', 'O2', 'GOC', 'PAR', 'PHY', 'so']

for varname in varlist:

    print(varname)
    
    data = xr.open_dataset('%s/hovmoller_data_0-%d_%s.nc' %(dirout, zmax, varname))
    oope = data[varname].to_masked_array()
    
    clim, anoms = ts.get_monthly_clim(oope)

    dataout = data.copy()
    dataout[varname].values = anoms

    fileout = '%s/hovmoller_data_0-%d_%s_anoms.nc' %(dirout, zmax, varname)
    dataout.to_netcdf(fileout)
    dataout.attrs['file'] = os.path.realpath(__file__)
    dataout.attrs['date'] = str(datetime.today())


