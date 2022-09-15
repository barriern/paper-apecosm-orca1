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
#     display_name: Python [conda env:nbarrier2]
#     language: python
#     name: conda-env-nbarrier2-py
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

ytt = [1982, 1997, 2015]
mtt = [10, 11, 12]
# -

oope0 = xr.open_dataset('/home1/scratch/nbarrier/test_pacific_ORCA1_JRAC02_CORMSK_CYC3_FINAL_OOPE_init_Y1958D000.nc')
oope0 = oope0['OOPE_init'].squeeze()
oope0 = float(oope0)


def process_time_series(varname):
    oope = xr.open_mfdataset('/home1/scratch/nbarrier/test_pacific_*%s*.nc' %varname)[varname].squeeze()
    oope

    dates = oope['time.year'] * 100 + oope['time.month']

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
    
    return output, dates


# +
varnames = [
    'madv_trend',
    'zadv_trend',
    'mdiff_trend',
    'zdiff_trend',
    'predationTrend',
    'growthTrend',
]

output = {}

for v in ['OOPE'] + varnames:
    temp, dates = process_time_series(v)
    output[v] = temp
output.keys()
# -

plt.figure()
oope_clim = output['OOPE']['clim']
oope_anoms = output['OOPE']['anom']
oope_tot = oope_anoms + oope_clim
plt.plot(oope_tot - oope_tot.mean())
plt.plot(oope_anoms)
plt.plot(oope_clim - oope_clim.mean())

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
list_dates = [y * 100 + m for y in ytt for m in mtt]
list_dates

test_dates = [d in list_dates for d in dates[1:].values]
idates = np.nonzero(test_dates)[0]
dates[1:][idates]

# +
series = {}
for v in varnames:
    
    temp = output[v]['anom']
    temp = np.cumsum(temp)
    temp = 0.5 * (temp[1:] + temp[:-1])
    series[v] = temp

plt.figure(figsize=(12, 8))
cpt = 1
for v in varnames:
    plt.subplot(3, 2, cpt)
    plt.plot(series[v], marker='.')
    plt.plot(idates, series[v][idates], marker='o', linestyle='none')
    plt.title(v)
    cpt += 1
# -

compo_oope_anom = output['OOPE']['anom'][1:][idates].mean()
compo_oope_anom

# +
toto = 0
for t in series:
    print(t)
    toto += series[t]
    
plt.figure(figsize=(12, 8))
plt.plot(output['OOPE']['anom'][1:])
plt.plot(toto)
plt.plot(idates, output['OOPE']['anom'][1:][idates], marker='.', linestyle='none')
plt.plot(idates, toto[idates], marker='.', linestyle='none')

toto[idates].mean(), output['OOPE']['anom'][1:][idates].mean()
# -


