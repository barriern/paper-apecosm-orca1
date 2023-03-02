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

grid = 'ptrc_T'
varlist = ['PHY2', 'ZOO', 'ZOO2', 'GOC']
# -

dirmd = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/'
mesh = xr.open_dataset('%s/mesh_mask_eORCA1_v2.2.nc' %dirmd).isel(t=0)
mesh

mesh = mesh.rename({'z': 'olevel'})
mesh

lon = mesh['glamt']
lat = mesh['gphit']
domain = (abs(lat) <= latmax)
domain = domain & ((lon >= lonwest) | (lon <= loneast))
domain

z1d = mesh['gdept_1d']
depth = (z1d <= zmax)
depth

filelist = glob('%s/JRA_CO2/nico*%s*nc' %(dirmd, grid))
filelist.sort()
filelist

for f in filelist[:]:
    
    print('++++++++++++++++++++++++ Processing file ', f)

    basename = os.path.basename(f)
    outfile = os.path.join(os.getenv('SCRATCH'), 'pacific_' + basename)
    outfile

    data = xr.open_dataset(f, decode_times=False)[varlist]
    data

    dataout = data.where((domain == True) & (depth == True), drop=True)
    dataout

    dataout.attrs['loneast'] = loneast
    dataout.attrs['lonwest'] = lonwest
    dataout.attrs['zmax'] = zmax
    
    plk = dataout['ZOO2'] + dataout['PHY2'] + dataout['ZOO'] + dataout['GOC']
    plk.name = 'PLK'
    
    encoding = {}
    encoding['PLK'] = {'zlib': True, "complevel": 9}
    
    plk.to_netcdf(outfile, encoding=encoding, unlimited_dims='time_counter')


