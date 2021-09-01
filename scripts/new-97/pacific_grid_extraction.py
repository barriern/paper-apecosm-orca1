# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Python [conda env:nbarrier] *
#     language: python
#     name: conda-env-nbarrier-py
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt
import os
from dask.diagnostics import ProgressBar

latmax = 42
lonwest = 110
loneast = -60
zmax = 250
# -

dirmd = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/'

mesh = xr.open_dataset('%s/mesh_mask_eORCA1_v2.2.nc' %dirmd).isel(t=0)
mesh

lon = mesh['glamt']
lon

lat = mesh['gphit']
lat

domain = (abs(lat) <= latmax)
domain = domain & ((lon >= lonwest) | (lon <= loneast))
domain

meshout = mesh.where(domain == True, drop=True)
meshout

z1d = mesh['gdept_1d']
z1d

depth = z1d <= zmax
depth

meshout = meshout.where(depth == True, drop=True)
meshout

meshout.isel(z=0)['tmask'].plot()

meshout.attrs['loneast'] = loneast
meshout.attrs['lonwest'] = lonwest
meshout.attrs['zmax'] = zmax

outdir = os.getenv('SCRATCH')
outfile = '%s/pacific_mesh_mask.nc' %outdir
outfile

with ProgressBar():
    meshout.to_netcdf(outfile)


