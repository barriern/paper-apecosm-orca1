# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python [conda env:nbarrier]
#     language: python
#     name: conda-env-nbarrier-py
# ---

import xarray as xr
from glob import glob
import numpy as np

# +
for varname in ['madv_trend', 'zadv_trend', 'mdiff_trend', 'zdiff_trend', 'growthTrend', 'predationTrend']:

    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ', varname)
    clim = xr.open_dataset('data/pacific_clim_%s.nc' %varname)[varname]
    clim

    filelist = []
    for y in [1982, 1997, 2015]:
        subfile = []
        pattern = '/home1/scratch/nbarrier/pacific_ORCA1*%s*%d*.nc' %(varname, y)
        temp = glob(pattern)[0]
        subfile.append(temp)
        pattern = '/home1/scratch/nbarrier/pacific_ORCA1*%s*%d*.nc' %(varname, y + 1)
        temp = glob(pattern)[0]
        subfile.append(temp)
        pattern = '/home1/scratch/nbarrier/pacific_ORCA1*%s*%d*.nc' %(varname, y - 1)
        temp = glob(pattern)[0]
        subfile.append(temp)
        subfile.sort()
        filelist.append(subfile)
    filelist.sort()
    filelist

    output = 0
    for f in filelist:
        print(f)
        data = xr.open_mfdataset(f)[varname]
        anom = data.groupby('time.month') - clim
        print(anom['time'].values)
        output += anom.values
    output /= len(filelist)
    output.shape

#     Output is the mean value of trend anomalies over the 3 years (36 months) arount the El Nino. Dimensions are `(nmonths=36, ny, nx, nl)`

#     Now we will try to extract a seasonal kind time-series.

    time = np.arange(output.shape[0])
    time

    season = time // 3
    season

    dsout = xr.Dataset()
    dsout[varname] = (['time', 'y', 'x', 'l'], output)
    dsout['season'] = (['time'], season)
    dsout

    output = dsout.groupby(dsout['season']).mean(dim='time')
    output

    foutname = 'data/seasonal_composite_trends_%s.nc' %varname
    foutname

    output.to_netcdf(foutname)
    
    del(clim)
    del(output)
# -


