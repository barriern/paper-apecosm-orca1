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
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid, ImageGrid

species = 'SKJ'
# -

# ## Processing Apecosm outputs

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc').rename({'wpred': 'w'})
length = const['length'].values * 100
const = const.assign(w=length)

data = xr.open_dataset('apecosm_time_series.nc')
data = data.assign(w=length)
data

apeeast = (data['tseast'] * const['weight_step']).sel(w=slice(25, 70)).sum(dim='w')
apeeast

apewest = (data['tswest'] * const['weight_step']).sel(w=slice(25, 70)).sum(dim='w')
apewest

# ## Processing Sardara outputs

data = xr.open_dataset('time_series_gear_%s_species_PS_latmax_2.nc' %species)
data

dates = data['time'].values
years = dates // 100
years

iok = np.nonzero((years >= 1958) & (years <= 2018))[0]
iok

data = data.isel(time=iok)
data

sareast = data['tseast'].values
sarwest = data['tswest'].values

dates = data['time'].values
years = dates // 100
months = dates - 100 * years
labels = np.array(['%.4d-%.2d' %(y, m) for y, m in zip(years, months)])


#
# ## Plotting outputs

def compute_anoms(data):
    
    index = np.arange(12)
    nyears = len(data) // 12
    clim = np.zeros((12))
    for i in range(nyears):
        clim += data[index]
        index += 12
    clim /= nyears
    
    index = np.arange(12)
    nyears = len(data) // 12
    anom = np.zeros(data.shape)
    for i in range(nyears):
        anom[index] = data[index] - clim
        index += 12
    
    return anom


iok = np.nonzero(years >= 1990)[0]
if(True): 
    sareast = compute_anoms(sareast[iok])
    sarwest = compute_anoms(sarwest[iok])
    apeeast = compute_anoms(apeeast.values[iok])
    apewest = compute_anoms(apewest.values[iok])
else:
    sareast = sareast[iok]
    sarwest = sarwest[iok]
    apeeast = apeeast.values[iok]
    apewest = apewest.values[iok]

# +
fig = plt.figure(facecolor='white', figsize=(12, 12))
plt.rcParams['font.size'] = 15

#axgr = ImageGrid(fig, 111,  nrows_ncols=(2, 1), 
#                 label_mode='L', aspect=False, share_all=False, axes_pad=[1, 0.5])

ax = plt.subplot(211)
plt.plot(sareast)
ax2 = ax.twinx()
plt.plot(apeeast, color='firebrick')
ax2.set_ylabel('OOPE (J)', color='firebrick')
ax.set_ylabel('%s catch' %species)
plt.title('East')
ax2.spines['right'].set_color('FireBrick')
plt.setp(ax2.get_yticklabels(), color='firebrick')
ax2.tick_params(color='firebrick')
ax2.yaxis.get_offset_text().set(color='firebrick')
plt.setp(ax.get_xticklabels(), visible=False)
ax.grid(True)

ax = plt.subplot(212, sharex=ax)
plt.plot(sarwest)
ax2 = ax.twinx()
plt.plot(apewest, color='firebrick')
ax2.set_label('OOPE')
plt.title('West')
ax2.spines['right'].set_color('FireBrick')
plt.setp(ax2.get_yticklabels(), color='firebrick')
ax2.tick_params(color='firebrick')
ax2.yaxis.get_offset_text().set(color='firebrick')
time = np.arange(len(iok))
stride = 2 * 12
t = ax.set_xticks(time[::stride])
l = ax.set_xticklabels(labels[iok][::stride], rotation=45, ha='right')
ax.set_xlim(time.min(), time.max())
ax.grid(True)
ax2.set_ylabel('OOPE (J)', color='firebrick')
ax.set_ylabel('%s catch' %species)
plt.savefig('time_series_west_east_species_%s.png' %(species), bbox_inches='tight')
# -


