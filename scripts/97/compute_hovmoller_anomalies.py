# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Plotting equatorial Hovmoller

# ## Imports

# +
import xarray as xr
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = False

varname = 'thetao'
# -

# ## Domain extraction
#
# Here, the mesh file is masked when latitude is beyond 2N/2S 

latmax = 2
mesh = xr.open_dataset('extracted_mesh_file.nc').isel(z=0)
mesh

lat = mesh['gphit']
lat

mesh = mesh.where(abs(lat) <= latmax)
mesh

# Here, the surface is extracted and multiplied by the land-sea mask values.

surf = mesh['e1t'] * mesh['e2t'] * mesh['tmask']
surf

# Finally, the mean longitude is extracted. It will be used as a coordinate for the Hovmoller.

lonf = mesh['glamt'].mean(dim='y')
lonf.plot()

# The longitudes are now modified to match the Pacific system (0-360 instead of -180/180)

lonf = (lonf + 360) % 360

# ## Computation of Hovmoller anomalies

anoms = xr.open_dataset('anom_%s.nc' %varname)[varname]
anoms

# Presumably, no need to mask anomalies since surface is masked already.

output = (anoms * surf).sum(dim='y') / (surf.sum(dim='y'))
output

output['x'] = lonf.values
output

output = output.sortby(output['x'])
output

output.plot()

output.to_netcdf('hov_anomalies_%s.nc' %varname)


