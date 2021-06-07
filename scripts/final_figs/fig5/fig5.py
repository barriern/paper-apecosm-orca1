import cartopy.mpl.gridliner as gridliner
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os.path
import datetime
import sys
import scipy.signal as sig
from mpl_toolkits.axes_grid1.axes_grid import ImageGrid

import string
letters = list(string.ascii_lowercase)

    
def get_data(varname):
    
    data = xr.open_dataset('data/equatorial_covariance_%s.nc' %varname)
    lon = data['lon']
    lon = (lon + 360) % 360
    data['lon'] = lon

    data = data.sortby(data['lon'])
    lon = data['lon'].values
    lags = data['lags'].values
    cov = data['covariance'].values
    
    return lags, lon, cov

xsize = 280
xvar = 160
ysize = -50
ylet = 50

_DEGREE_SYMBOL = u'\u00B0'

def _east_west_formatted(longitude, num_format='g'):
    fmt_string = u'{longitude:{num_format}}{degree}{hemisphere}'
    if(longitude == 180):
        longitude = 0
    output = fmt_string.format(longitude=abs(longitude), num_format=num_format,
                             hemisphere=gridliner._lon_hemisphere(longitude),
                             degree=_DEGREE_SYMBOL)
    return output

import re
from glob import glob

pattern = 'data/equatorial_covariance_(.*).nc'
regexp = re.compile(pattern)
filelist = glob('data/*equatorial*nc')
filelist.sort()
varlist = []
for f in filelist:
    varname = regexp.match(f).groups()[0]
    varlist.append(varname)

constant = xr.open_dataset('../../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
constant = constant.isel(wpred=[14, 45, 80])
wstep = constant['weight_step'].values
length = constant['length'].values * 100

sizes = ['%.dcm' %l for l in length]
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)


############################## extract data to plot
lags, lon, oope = get_data('OOPE')
oope = oope * wstep[:, np.newaxis, np.newaxis]

lags, lon, repfonct = get_data('repfonct_day')  # w, lon, time
lags, lon, mort_day = get_data('mort_day')  # w, lon, time
lags, lon, gamma1 = get_data('gamma1')  # w, lon, time
lags, lon, zadv_trend = get_data('zadv_trend')  # w, lon, time
lags, lon, madv_trend = get_data('madv_trend')  # w, lon, time
zadv_trend = zadv_trend * wstep[:, np.newaxis, np.newaxis]
madv_trend = madv_trend * wstep[:, np.newaxis, np.newaxis]

varlist = [oope, repfonct, mort_day, gamma1, zadv_trend, madv_trend]
varlist = [oope, repfonct, mort_day, gamma1]
varnames = ['oope', 'repfonct', 'mortday', 'gamma1', 'zadv', 'madv']

nvar = len(varlist)
    
fig = plt.figure(figsize=(14, 8))
axgr = ImageGrid(fig, 111, nrows_ncols=(nvar, 3), ngrids=None, direction='row', axes_pad=(0.8, 0.2), share_all=False, aspect=True, label_mode='L', cbar_mode='each', cbar_location='right', cbar_pad='5%', cbar_size='5%', cbar_set_cax=True, axes_class=None)
cbar_axes = axgr.cbar_axes
axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

cpt = 0

for c in range(nvar):  # loop over the variables

    varname = varnames[c]
    cov = varlist[c]
    print('Plotting variable ', varname)
    print(lon)

    mask = (lon >= 130) & (lon <= 300)
    ilon = np.nonzero(mask == False)[0]
    ilags = np.nonzero(np.abs(lags) <= 3*12)[0]

    ntime = (61 * 12 - 1)

    for i in range(3):
        covtemp = cov[i]
        covtemp[ilon, :] = np.ma.masked
        #covtemp[:, ilags] = np.ma.masked
        covtemp = np.ma.masked_where(covtemp == 0, covtemp)
        iok = np.nonzero(covtemp.mask == False)
        covtemp /= ntime

        ax = axout[cpt]

        perc = np.percentile(np.ravel(np.abs(covtemp[iok])), 99.5)
        levels = np.linspace(-perc, perc, 11)
        cs = ax.pcolormesh(lon, lags, covtemp.T, shading='auto')
        cl = ax.contour(lon, lags, covtemp.T, levels=levels, colors='k', linewidths=0.5)
        cb = axgr.cbar_axes[cpt].colorbar(cs)
        ##axgr.cbar_axes[cpt].xaxis.set_ticks_position('top')
        cs.set_clim(-perc, perc)
        ax.set_xlim(130, 300)
        ax.grid(True)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Lags (months)')
        #if(varname == 'OOPE'):
        #    cb.set_label('J/m2')
        #ax.set_title('%d cm' %length[i])

        labels = ['150', '180', '-150', '-120', '-90', '-60']
        xticks = np.array([float(l) for l in labels])
        xticks[xticks < 0] += 360
        labels = [_east_west_formatted(float(l)) for l in labels]

        ax.text(xsize, ysize, sizes[i], bbox=dicttext, ha='center', va='center')
        ax.text(xvar, ysize, varname, bbox=dicttext, ha='center', va='center')
        ax.text(xsize, ylet, letters[cpt] + ')', bbox=dicttext, ha='center', va='center')

        ax.set_xticklabels(labels)
        ax.set_xticks(xticks)
        cpt += 1

    plt.savefig('fig5', bbox_inches='tight')

#mask = (lon >= 130) & (lon <= 300)
#ilon = np.nonzero(mask == True)[0]
#
#plt.figure()
#plt.subplot(211)
#iii = ilon[0]
#plt.plot(zadv_trend[0, iii, :])
#plt.title('zadv-trend')
#
#plt.subplot(212)
#iii = ilon[0]
#plt.plot(izadv_trend[0, iii, :])
#plt.title('cumsum zadv-trend')
#
#plt.savefig('ts')
#
