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

# +
import xarray as xr
import matplotlib.pyplot as plt
import os
from dask.diagnostics import ProgressBar
from glob import glob

latmax = 42
lonwest = 110
loneast = -60
zmax = 250

ltoe = [3, 20, 90]

grid = 'OOPE_init'
# -

dirmd = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/'
mesh = xr.open_dataset('%s/mesh_mask_eORCA1_v2.2.nc' %dirmd).isel(t=0)
mesh

mesh = mesh.rename({'z': 'depth'})
mesh

lon = mesh['glamt']
lat = mesh['gphit']
domain = (abs(lat) <= latmax)
domain = domain & ((lon >= lonwest) | (lon <= loneast))
domain

z1d = mesh['gdept_1d']
depth = (z1d <= zmax)
depth

dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output'
pattern = '%s/*%s*nc' %(dirin, grid)
print(pattern)
filelist = glob(pattern)
filelist.sort()
filelist = filelist[:]
filelist[:5]

pattern = '%s/*%s*nc' %(dirin, 'Constant')
print(pattern)
constfile = glob(pattern)[0]
constfile

const = xr.open_dataset(constfile)
const

length = const['length'] * 100
length = length.rename({'wpred': 'w'})
length

data = xr.open_dataset(filelist[0])
data['w'] = length
data

if(grid == 'FORAGE'):
    dataout = data.where((domain == True) & (depth == True), drop=True)
else:
    dataout = data.where((domain == True), drop=True)
dataout

for f in filelist:
    
    print('++++++++++++++++++++++++ Processing file ', f)

    basename = os.path.basename(f)
    outfile = os.path.join(os.getenv('SCRATCH'), '3size_pacific_' + basename)
    outfile
    if os.path.isfile(outfile):
        print('File %s exists. Skipped' %outfile)
        continue

    data = xr.open_dataset(f, decode_times=False).isel(community=0)
    if(grid != 'FORAGE'):
        data['w'] = length
        data = data.sel(w=ltoe, method='nearest')
        data

    if(grid == 'FORAGE'):
        dataout = data.where((domain == True) & (depth == True), drop=True)
    else:
        dataout = data.where((domain == True), drop=True)
    dataout

    dataout.attrs['loneast'] = loneast
    dataout.attrs['lonwest'] = lonwest
    dataout.attrs['zmax'] = zmax

    print('writting ', outfile)
    
    varname = [v for v in dataout.variables if v != 'time'][0]
    print(varname)
    encode = {varname: {'zlib': True, "complevel": 9}}
    
    dataout.to_netcdf(outfile, encoding=encode, unlimited_dims='time')


