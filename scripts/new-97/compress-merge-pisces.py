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
#     display_name: Python [conda env:nbarrier2]
#     language: python
#     name: conda-env-nbarrier2-py
# ---

# +
import xarray as xr
from glob import glob
from dask.diagnostics import ProgressBar
import os

grid = 'ptrc_T'

varlist = {}
varlist['grid_T'] = ['thetao', 'e3t']
varlist['ptrc_T'] = ['ZOO2','ZOO','GOC', 'PHY2']
varlist['add_T'] = ['NCHL', 'DCHL']
varlist['speed_V'] = ['vo']
varlist['speed_U'] = ['uo']

pattern = '/home1/scratch/nbarrier/pacific_nico_GCB-eORCA1-JRA14-CO2ANTH_*%s.nc' %grid
pattern

filelist = glob(pattern)
filelist.sort()
filelist

comp = dict(zlib=True, complevel=9)
encoding = {var: comp for var in varlist[grid]}
encoding

print(encoding)

for f in filelist:
    
    print('++++++++++++++++ Processing file ', f)
    
    data = xr.open_mfdataset(filelist, decode_times=False)
    data

    if((grid == 'speed_U') | (grid == 'speed_V')):
        data = data.drop(['time_centered'])
    else:       
        data = data.drop(['time_centered', 'time_counter_bounds'])
    data

    data = xr.decode_cf(data)
    data['time_counter']

    data = data[varlist[grid]]
    data

    foutname = '/home1/scratch/nbarrier/compressed_%s' %(os.path.basename(f))
    foutname
    print(foutname)

    data.to_netcdf(foutname, encoding=encoding)
