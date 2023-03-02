import xarray as xr
import scipy.signal as sig
import numpy as np
import os.path
import scipy.stats as stats

# +
import sys
sys.path.append("../..")
from extract_nino import read_index

dirin = '.'
# -

daten, nino = read_index(filename='../../data/external/oni.data', keepnan=False)

data = xr.open_mfdataset("%s/anom_interpolated*nc" %dirin, combine='by_coords')
lat = data['lat'].values
lon = data['lon'].values
date = data['time.year'] * 100 + data['time.month']
date = date.values

dmax = np.min([date.max(), daten.max()])
dmin = np.max([date.min(), daten.min()])
dmax, dmin

inino = np.nonzero((daten <= dmax) & (daten >= dmin))
idata = np.nonzero((date <= dmax) & (date >= dmin))[0]

daten = daten[inino]
nino = nino[inino]
date = date[idata]
data = data.isel(time=idata)
data = data['chl_anom'].values

realtime = np.arange(0, len(date))

data = data.T  # time, lat, lon -> lon, lat, time
data.shape

nlon, nlat, ntime = data.shape

output = np.zeros((nlon, nlat))
count = np.zeros((nlon, nlat))

for i in range(nlon):
    for j in range(nlat):
        temp = data[i, j]
        iok = np.nonzero(np.isnan(temp) == False)[0]
        count[i, j] = len(iok) / len(temp) * 100
        if len(iok) > 0:
            slope, intercept, r_value, p_value, std_err = stats.linregress(realtime[iok], temp[iok])
            trend = slope * realtime + intercept
            temp -= slope
            output[i, j] = np.cov(temp[iok], nino[iok])[0, 1]

output = output.T  # lat, lon
count = count.T

dsout = xr.Dataset({'lon': lon, 'lat':lat})
dsout['cov'] = (['lat', 'lon'], output)
dsout['count'] = (['lat', 'lon'], count)
dsout.to_netcdf("%s/covariance_satellite_data.nc" %dirin)


