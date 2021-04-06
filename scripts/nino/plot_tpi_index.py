import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

raw = pd.read_csv("data/tpi.timeseries.ersstv5.data", index_col=0, skiprows=1, skipfooter=11, delim_whitespace=True, engine='python', na_values=-99, header=None)
years = raw.index.values
month = np.arange(12) + 1
raw = np.ravel(raw.values)

fil = pd.read_csv("data/tpi.timeseries.ersstv5.filt.data", index_col=0, skiprows=1, skipfooter=11, delim_whitespace=True, engine='python', na_values=-99, header=None)
years = fil.index.values
month = np.arange(12) + 1
fil = np.ravel(fil.values)

month, years = np.meshgrid(month, years)
date = years * 100 + month
date = np.ravel(date)

iok = np.nonzero((date >= 195801) & (date <= 201812))[0]
date = date[iok]
fil = fil[iok]
raw = raw[iok]

plt.figure()
plt.plot(raw, label='raw')
plt.plot(fil, label='fil')
plt.savefig('tpi.png')

dsout = xr.Dataset()
dsout['tpi_filt'] = (['time'], fil)
dsout['tpi_raw'] = (['time'], raw)
dsout['time'] = (['time'], date)
dsout.to_netcdf('filt_tpi.nc')
