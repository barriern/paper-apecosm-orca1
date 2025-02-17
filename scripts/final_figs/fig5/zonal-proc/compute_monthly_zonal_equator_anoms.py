import xarray as xr
import numpy as np
import apecosm.ts as ts
import os.path 
from datetime import datetime

dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'
dirmesh = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST'

prefix = 'final-runs'

#for varname in ['mort_day']:
#for varname in ['gamma1', 'mort_day', 'mdiff_trend', 'madv_trend', 'zdiff_trend', 'zadv_trend', 'repfonct_day', 'OOPE']:
for varname in ['adv_trend']:

    filein = '%s/%s_%s_meridional_mean.nc' %(dirout, prefix, varname)
    try:
        data = xr.open_dataset(filein)
        oope = data[varname].values

        clim, anoms = ts.get_monthly_clim(oope)

        dataout = data.copy()
        dataout[varname].values = anoms

        fileout = '%s/%s_%s_meridional_mean_anoms.nc' %(dirout, prefix, varname)
        dataout.to_netcdf(fileout)
        dataout.attrs['file'] = os.path.realpath(__file__)
        dataout.attrs['date'] = str(datetime.today())
    except:
        pass
