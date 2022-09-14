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
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal

import sys
sys.path.append('../../nino')
from extract_nino import read_index
# -

ymin = 1997
mmin = 9
datemin = '%.4d-%.2d-01' %(ymin, mmin)
datemin

ymax = 2018
mmax = 12
datemax = '%.4d-%.2d-01' %(ymax, mmax)
datemax

sldates = slice(datemin, datemax)

# ## Processing Chl

chl = xr.open_dataset('data/obs_equatorial_mean.nc')
chl = chl['chl']
chl

# ## Reading equatorial zonal velocity

uofile = 'data/satellite_nino_34_uo.nc'
uo = xr.open_dataset(uofile)
uo = uo['uo']
uo = uo.sel(time=slice(datemin, datemax))

# ## Computing NIno34 SSH

ssh = xr.open_dataset('data/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M-m_1641569308380.nc')
ssh = ssh['sla']
ssh

lat = ssh['latitude']
lon = ssh['longitude']

latmin = -5
latmax = 5
lonmax = -120 + 360
lonmin = -170 + 360
lonmin, lonmax
test = (lat <= latmax) & (lat >= latmin)
test = test & (lon<=lonmax) & (lon>=lonmin)
sshnino = ssh.where(test)
sshnino.isel(time=0).plot()

sshnino = sshnino.mean(dim=['longitude', 'latitude'])
sshnino

# ## Process Nino index

dnino, nino = read_index('../../data/index/oni.data')
iok = np.nonzero((dnino >= (ymin * 100 + mmin)) & (dnino <= (ymax * 100 +mmax)))[0]
dnino = dnino[iok]
nino = nino[iok]
dnino.shape

# ## Computation of correlations

chlclim = chl.groupby('time.month').mean(dim='time')
chlanom = chl.groupby('time.month') - chlclim
chlanom = chlanom.sel(time=sldate)
chlanom

sshninoclim = sshnino.groupby('time.month').mean(dim='time')
sshninoanom = sshnino.groupby('time.month') - sshninoclim
sshninoanom.values = scipy.signal.detrend(sshninoanom.values)
sshninoanom = sshninoanom.sel(time=sldate)
sshninoanom

# +
plt.figure()
plt.subplot(311)
plt.plot(nino)
plt.subplot(312)
plt.plot(-chlanom)

plt.subplot(313)
plt.plot(sshninoanom)
# -

print(np.corrcoef(sshninoanom, chlanom)[0, 1])


