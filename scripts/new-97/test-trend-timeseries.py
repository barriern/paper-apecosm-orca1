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
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

oope0 = xr.open_dataset('data/test/test_pacific_ORCA1_JRAC02_CORMSK_CYC3_FINAL_OOPE_init_Y1958D000.nc')
oope0 = oope0['OOPE_init'].squeeze()
oope0 = float(oope0)


def process_time_series(varname):
    oope = xr.open_mfdataset('data/test/*%s*.nc' %varname)[varname].squeeze()
    oope
    datestart = '1971-01-01'
    dateend = '2000-12-31'

    oope_clim = oope.groupby('time.month').mean(dim='time')
    oope_clim
    
    oope_anoms = oope.groupby('time.month') - oope_clim
    oope_anoms = oope_anoms.values
    oope_clim = oope_clim.values
    
    nyears = len(oope_anoms) // (len(oope_clim))

    oope_clim = np.ravel(np.tile(oope_clim, (nyears, 1)))
    
    output = {}
    output['clim'] = oope_clim
    output['anom'] = oope_anoms
    
    return output


# +
varnames = [
    'growthTrend',
    'madv_trend',
    'mdiff_trend',
    'predationTrend',
    'zadv_trend',
    'zdiff_trend'
]

output = {}

for v in ['OOPE'] + varnames:
    temp = process_time_series(v)
    output[v] = temp
output.keys()
# -

plt.figure()
oope_clim = output['OOPE']['clim']
oope_anoms = output['OOPE']['anom']
oope_tot = oope_anoms + oope_clim
plt.plot(oope_tot)
plt.plot(oope_clim + oope_anoms)

test = 0
for v in varnames[:]:
    print(v)
    test += output[v]['clim'] + output[v]['anom']
test.shape

# +
trend = oope0 + np.cumsum(test)
trend = 0.5 * (trend[1:] + trend[:-1])
print(trend.shape, oope_tot[1:].shape)
print(np.corrcoef(trend, oope_tot[1:]))

plt.figure(figsize=(12, 8))
plt.plot(trend)
plt.plot(oope_tot[1:])
# -

test = 0
for v in varnames[:]:
    print(v)
    print(np.mean(output[v]['anom']))
    test += output[v]['anom']
test.shape

# +
test = test
trend = np.cumsum(test)
trend = 0.5 * (trend[1:] + trend[:-1])

oope_to_plot = output['OOPE']['anom'][1:]
print(trend.shape, oope_to_plot.shape)

print(trend.shape, oope_tot[1:].shape)
print(np.corrcoef(trend, oope_tot[1:]))

plt.figure(figsize=(12, 8))
plt.plot(output['OOPE']['anom'])
plt.plot(trend)
# -


