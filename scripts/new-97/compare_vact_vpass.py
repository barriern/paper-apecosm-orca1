# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import string
from cartopy.mpl.ticker import LongitudeFormatter

zmax = 50
formatter0 = LongitudeFormatter(dateline_direction_label=True)

latmax = 5

letters = list(string.ascii_lowercase)

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
# -

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100
const

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(y=61)
mesh = mesh.rename({'z': 'olevel'})
mesh['olevel'] = mesh['gdept_1d'].isel(x=0)
mesh = mesh.where(mesh['olevel'] <= zmax)

lon0 = mesh['glamt'].isel(olevel=0)
lon0 = (lon0 + 360) % 360
lon0

boolean = ((lon0 >= 150) & (lon0 <= 270))
lon0 = lon0.where(boolean, drop=True)

lengths = [3, 20, 90]
lengths

data = xr.open_mfdataset('data/nino_equatorial_composites_*.nc', compat='override').where(boolean, drop=True)
data['l'] = const['l']
data

# ## Plotting the hovmoller

x = lon0.values
x

months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
tlabels = ['%s-%d' %(m, y) for y in [0, 1] for m in months]
tlabels

# +
max_vpass = abs(data['v_passive']).max(dim=['time', 'x'])
max_vpass

max_upass = abs(data['u_passive']).max(dim=['time', 'x'])
max_upass

u_pass = xr.concat([max_upass, max_vpass], dim='d')
u_pass = u_pass.max(dim='d')
u_pass.plot()

# +
max_vact = abs(data['v_active']).max(dim=['time', 'x'])
max_vact

max_uact = abs(data['u_active']).max(dim=['time', 'x'])
max_uact

u_act = xr.concat([max_uact, max_vact], dim='d')
u_act = u_act.max(dim='d')
u_act.plot()
# -

factor = (u_pass / u_act)
ax = plt.gca()
factor.plot(ax=ax)
ax.set_yscale('log')
factor.sel(l=lengths, method='nearest').values

i = 2
plt.figure(figsize=(18,18))
plt.subplot(221)
data['u_active'].sel(l=lengths[i], method='nearest').plot()
plt.subplot(222)
data['v_active'].sel(l=lengths[i], method='nearest').plot()
plt.subplot(223)
data['u_passive'].sel(l=lengths[i], method='nearest').plot()
plt.subplot(224)
data['v_passive'].sel(l=lengths[i], method='nearest').plot()


