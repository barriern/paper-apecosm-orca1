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

# # Plotting time-series
#
# ## Imports

# +
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

varname = 'adv_trend'
# -

# ## Loading constant fields

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const

# ## Plotting time-series

# +
dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/processed_pacific'

data = xr.open_mfdataset('%s/*_timeseries.nc' %(dirin))
data
# -

data = data.rename({'w': 'length'})
data

data.coords['length'] = const['length'].values * 100
data['length']

data = data.sel(length=[3, 20, 90], method='nearest').sel(time=slice('1997-01-01', '1998-12-31'))
data

if(varname == 'adv_trend'):
    print('Computing the sum of zonal and meriodonal advective trends + time-integration')
    var_central = data['madv_trend_central'].values + data['zadv_trend_central'].values
    var_central = np.cumsum(var_central, axis=0) * 24 * 60 * 60 * 30
    
    var_west = data['madv_trend_west'].values + data['zadv_trend_west'].values
    var_west = np.cumsum(var_west, axis=0) * 24 * 60 * 60 * 30
else:
    var_central = data['%s_central' %varname].values
    var_west = data['%s_west' %varname].values

const = const.rename({'wpred': 'l'})
const['l'] = const['length'].values * 100
const = const.sel(l=[3, 20, 90], method='nearest')
const

date = data['time'].values

years = np.array([d.year for d in date])
months = np.array([d.month for d in date])
datestr = ['%.4d-%.2d' %(y, m) for y, m in zip(years, months)]
time = np.arange(len(years))

if (varname == 'OOPE') | (varname == 'adv_trend'):
    print('Convert in J/m2')
    var_central *= const['weight_step'].values
    var_west *= const['weight_step'].values

length = const['length'].values * 100
length

# +
fig = plt.figure(figsize=(8, 14))

plt.rcParams['font.size'] = 14

stride = 2

# size=3cm
i = 1
ax = plt.subplot(3, 1, i)
plt.plot(time, var_central[:, i - 1], label='central')
plt.plot(time, var_west[:, i - 1], label='west')
ax.set_title('%.f cm' %length[i - 1])
plt.legend(loc = 0)
ax.set_ylabel(varname)
ax.get_xaxis().set_visible(False)  # removes xlabels
ticks = ax.set_xticks(time[::stride])
ylim = ax.get_ylim()
yyy = np.abs(ylim).max()
ax.set_ylim(-yyy, yyy)
ax.set_xlim(0, time.max())

# size=3cm
i = 2
ax = plt.subplot(3, 1, i)
plt.plot(time, var_central[:, i - 1], label='central')
plt.plot(time, var_west[:, i - 1], label='west')
ax.set_title('%.f cm' %length[i - 1])
ax.set_ylabel(varname)
ticks = ax.set_xticks(time[::stride])
ax.get_xaxis().set_visible(False)  # removes xlabels
ylim = ax.get_ylim()
yyy = np.abs(ylim).max()
ax.set_ylim(-yyy, yyy)
ax.set_xlim(0, time.max())

i = 3
ax = plt.subplot(3, 1, i)
plt.plot(time, var_central[:, i - 1], label='central')
plt.plot(time, var_west[:, i - 1], label='west')
ax.set_title('%.f cm' %length[i - 1])
ax.set_ylabel(varname)
ticks = ax.set_xticks(time[::stride])
lab = ax.set_xticklabels(datestr[::stride], ha='right', rotation=45)
ylim = ax.get_ylim()
yyy = np.abs(ylim).max()
ax.set_ylim(-yyy, yyy)
ax.set_xlim(0, time.max())

plt.savefig('time-series_%s.pdf' %varname, bbox_inches='tight')
# -


