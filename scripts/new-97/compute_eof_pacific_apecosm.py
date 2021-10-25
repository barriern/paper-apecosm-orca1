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

# # EOF calculation
#
# In this section, EOF calculation is performed on detrended files.

# ## Imports

# +
import xarray as xr
from eofs.xarray import Eof
import os
from glob import glob
import numpy as np
from dask.diagnostics import ProgressBar

latmax = 10
lonmin = 150
lonmax = -100
# -

# ## Computes weights

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0)
mesh

lon = mesh['glamt']
lon

lat = mesh['gphit']
lat

weights = mesh['e1t'] * mesh['e2t'] * mesh['tmask'] * mesh['tmaskutil']
weights

weights.plot()

domain = abs(lat) <= latmax
domain = domain & ((lon >= lonmin) | (lon <= lonmax))
weights = weights.where(domain)
weights.plot()

# Now the weights are normalized:

weights = weights / (weights.sum())
weights.plot()

weights = np.sqrt(weights)
weights.plot()

# ## Computation of the EOFs

# First, the data are loaded into memory

dirin = os.getenv('SCRATCH')
data = xr.open_dataset(os.path.join(dirin, 'detrended_pacific_OOPE_anomalies.nc'))
data = data['OOPE']
data

ny, nx, nw, ntime = data.shape
ny, nx, nw, ntime

# Init output values:

eofmaps = np.zeros((nw, 2, ny, nx), dtype=float)
eofcovmaps = np.zeros((nw, 2, ny, nx), dtype=float)
eofpcs = np.zeros((nw, ntime, 2), dtype=float)
eofvar = np.zeros((nw, 2))
eofmaps.shape

for l in range(0, nw):
    print('Processing w = ', l)
    temp = data.isel(w=l)
    temp

    solver = Eof(temp.T, weights=weights.T)

    # Now, computation of EOF maps as covariance maps

    covmaps = solver.eofsAsCovariance(neofs=2).to_masked_array()  # eof, x, y,
    eofcovmaps[l] = np.transpose(covmaps, (0, 2, 1))

    # Computation of EOF maps as non-dimensionalized EOFs

    maps = solver.eofs(eofscaling=0, neofs=2)
    eofmaps[l] = np.transpose(covmaps, (0, 2, 1))
    
    eofvar[l, :] = solver.varianceFraction(neigs=2) * 100

    # Computation of PCs
    
    pcs = solver.pcs(pcscaling=1, npcs=2)  # ntime, neofs
    eofpcs[l] = pcs

dsout = xr.Dataset()
dsout['eofcovmaps'] = (['w', 'eof', 'y', 'x'], eofcovmaps)
dsout['eofmaps'] = (['w', 'eof', 'y', 'x'], eofmaps)
dsout['eofpcs'] = (['w', 'time', 'eof'], eofpcs)
dsout['eofvar'] = (['w', 'eof'], eofvar)
dsout['time'] = data['time']
dsout.attrs['lonmin'] = lonmin
dsout.attrs['lonmax'] = lonmax
dsout.attrs['latmax'] = latmax
dsout

fileout = '%s/full_eof_pacific_OOPE_latmax_%d_lonmin_%d_lonmax_%d.nc' %(dirin, latmax, lonmin, lonmax)
fileout

with ProgressBar():
    dsout.to_netcdf(fileout)
