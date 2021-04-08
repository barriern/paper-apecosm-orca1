import cartopy.mpl.gridliner as gridliner
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os.path
import datetime
from cftime import utime
import sys
import scipy.signal as sig
from mpl_toolkits.axes_grid1.axes_grid import ImageGrid

_DEGREE_SYMBOL = u'\u00B0'

def _north_south_formatted(longitude, num_format='g'):
    fmt_string = u'{longitude:{num_format}}{degree}{hemisphere}'
    output = fmt_string.format(longitude=abs(longitude), num_format=num_format,
                             hemisphere=gridliner._lat_heimisphere(longitude),
                             degree=_DEGREE_SYMBOL)
    return output

import re
from glob import glob

lonmax = -150

pattern = 'data/meridional_covariance_%.f_(.*).nc' %(lonmax)
regexp = re.compile(pattern)
filelist = glob('data/*meridional*nc')
filelist.sort()
varlist = []
for f in filelist:
    varname = regexp.match(f).groups()[0]
    varlist.append(varname)

constant = xr.open_dataset('../../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
constant = constant.isel(wpred=[14, 45, 80])
wstep = constant['weight_step'].values
length = constant['length'].values * 100

for varname in varlist:

    print('Plotting variable ', varname)

    data = xr.open_dataset('data/meridional_covariance_%.f_%s.nc' %(lonmax, varname))
    lat = data['lat'].values
    lags = data['lags'].values
    cov = data['covariance']
    if varname == 'OOPE':
        cov = cov * wstep[:, np.newaxis, np.newaxis]
    cov = cov.values

    latmax = 40
    mask = (lat >= -latmax) & (lat <= latmax)
    ilon = np.nonzero(mask == False)

    fig = plt.figure(figsize=(14, 8))
    axgr = ImageGrid(fig, 111, nrows_ncols=(3, 1), ngrids=None, direction='row', axes_pad=(0.05, 0.3), share_all=False, aspect=True, label_mode='L', cbar_mode='each', cbar_location='right', cbar_pad='5%', cbar_size='5%', cbar_set_cax=True, axes_class=None)

    ntime = (61 * 12 - 1)

    for i, ax in enumerate(axgr):
        covtemp = cov[i]
        print(covtemp.shape)
        covtemp[ilon, :] = np.ma.masked
        #covtemp[:, ilags] = np.ma.masked
        covtemp = np.ma.masked_where(covtemp == 0, covtemp)
        iok = np.nonzero(covtemp.mask == False)
        covtemp /= ntime
        perc = np.percentile(np.ravel(np.abs(covtemp[iok])), 99.5)
        levels = np.linspace(-perc, perc, 11)
        cs = ax.pcolormesh(lat, lags, covtemp.T, shading='auto', cmap=plt.cm.RdBu_r)
        #cl = ax.contour(lat, lags, covtemp.T, levels=levels, colors='k', linewidths=0.5)
        cb = axgr.cbar_axes[i].colorbar(cs)
        cs.set_clim(-perc, perc)
        ax.set_xlim(-latmax, latmax)
        ax.grid(True)
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Lags (months)')
        #if(varname == 'OOPE'):
        #    cb.set_label('J/m2')
        ax.set_title('%d cm' %length[i])

        labels = ['150', '180', '-150', '-120', '-90', '-60']
        labels = np.arange(-latmax, latmax + 10, 10)
        xticks = np.array([float(l) for l in labels])
        labels = [_north_south_formatted(float(l)) for l in labels]

        ax.set_xticklabels(labels)
        ax.set_xticks(xticks)

    plt.suptitle(varname.replace('_', '-'))

    plt.savefig('meridional_covariance_%.f_%s.png' %(lonmax, varname), bbox_inches='tight')
