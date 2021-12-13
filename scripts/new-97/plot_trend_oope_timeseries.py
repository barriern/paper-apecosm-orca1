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
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import xarray as xr
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = False

domain = 'east'
l = 90

# ## Reading OOPE time-series

oope = xr.open_dataset('oope_timeseries.nc')
oope

tempoope = oope['oope_%s' %domain].sel(l=l, method='nearest')
tempoope

# ## Reading trends

trends = xr.open_mfdataset('*rend*timeseries.nc')
trends

temp = trends.sel(l=l, method='nearest')
temp

varnames = [v for v in trends.variables if domain in v]
varnames

for v in varnames[:]:
    print('v = ', v)
    if(v == varnames[0]):
        output = trends.sel(l=l, method='nearest')[v]
    else:
        output += trends.sel(l=l, method='nearest')[v]
output.name = 'trend'
output

# +
plt.figure(figsize=(12, 8))
plt.subplots_adjust(hspace=0.3)
plt.rcParams['font.size'] = 15
ax = plt.subplot(211)
(tempoope - tempoope.isel(time=0)).plot(label='OOPE')
(output - output.isel(time=0)).plot(label='trends')
ax.set_ylabel('J/m2')
ax.grid(True, linestyle='--')
plt.legend()

toto = trends.sel(l=l, method='nearest')['predationTrend_%s' %domain] + trends.sel(l=l, method='nearest')['growthTrend_%s' %domain] + trends.sel(l=l, method='nearest')['zadv_trend_%s' %domain] + trends.sel(l=l, method='nearest')['madv_trend_%s' %domain] + trends.sel(l=l, method='nearest')['zdiff_trend_%s' %domain] + trends.sel(l=l, method='nearest')['mdiff_trend_%s' %domain] 
ax2 = plt.subplot(212, sharex=ax)
trends.sel(l=l, method='nearest')['predationTrend_%s' %domain].plot(label='Pred.')
trends.sel(l=l, method='nearest')['growthTrend_%s' %domain].plot(label='Growth')
(trends.sel(l=l, method='nearest')['zadv_trend_%s' %domain] + trends.sel(l=l, method='nearest')['madv_trend_%s' %domain]).plot(label='Adv.')
(trends.sel(l=l, method='nearest')['zdiff_trend_%s' %domain] + trends.sel(l=l, method='nearest')['mdiff_trend_%s' %domain]).plot(label='Diff.')
toto.plot(label='Sum')
t = ax2.set_title('')
plt.legend(ncol=2)
ax2.grid(True, linestyle='--')
ax2.set_ylabel('J/m2')
plt.savefig('oope_trends_l_%d_domain_%s.png' %(l, domain), bbox_inches='tight')
# -

