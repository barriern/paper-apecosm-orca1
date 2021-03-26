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

pattern = 'data/final-runs_(.*)_meridional_mean_anoms.nc'
regexp = re.compile(pattern)
filelist = glob('data/*meridional_mean_anoms*nc')
filelist.sort()
varlist = []
for f in filelist:
    varname = regexp.match(f).groups()[0]
    varlist.append(varname)

varlist.remove('starvation')

constant = xr.open_dataset('/home/barrier/Work/apecosm/ORCA1/figures_orca1/plot_hov/ORCA1_JRA_CO2_CYC4_ConstantFields.nc')
constant = constant.isel(w=[14, 45, 80])
wstep = constant['weight_step'].values
length = constant['length'].values * 100

for varname in varlist:

    print('Plotting variable ', varname)

    data = xr.open_dataset('data/equatorial_covariance_%s.nc' %varname)
    lon = data['lon']
    lon = (lon + 360) % 360
    data['lon'] = lon

    data = data.sortby(data['lon'])
    lon = data['lon'].values
    lags = data['lags'].values
    cov = data['covariance']
    if varname == 'OOPE':
        cov = cov * wstep[:, np.newaxis, np.newaxis]
    cov = cov.values

    mask = (lon >= 130) & (lon <= 300)
    ilon = np.nonzero(mask == False)
    ilags = np.nonzero(lags < 0)[0]

    fig = plt.figure(figsize=(12, 8))
    axgr = ImageGrid(fig, 111, nrows_ncols=(3, 1), ngrids=None, direction='row', axes_pad=(0.05, 0.3), share_all=False, aspect=True, label_mode='L', cbar_mode='each', cbar_location='right', cbar_pad='5%', cbar_size='5%', cbar_set_cax=True, axes_class=None)

    ntime = (61 * 12 - 1)

    for i, ax in enumerate(axgr):
        covtemp = cov[i]
        covtemp[ilon, :] = np.ma.masked
        #covtemp[:, ilags] = np.ma.masked
        covtemp = np.ma.masked_where(covtemp == 0, covtemp)
        iok = np.nonzero(covtemp.mask == False)
        covtemp /= ntime

        perc = np.percentile(np.ravel(np.abs(covtemp[iok])), 99.9)
        levels = np.linspace(-perc, perc, 11)
        cs = ax.pcolormesh(lon, lags, covtemp.T, shading='auto')
        cl = ax.contour(lon, lags, covtemp.T, levels=levels, colors='k', linewidths=0.5)
        cb = axgr.cbar_axes[i].colorbar(cs)
        cs.set_clim(-perc, perc)
        ax.set_xlim(130, 300)
        ax.grid(True)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Lags (months)')
        if(varname == 'OOPE'):
            cb.set_label('J/m2')
        ax.set_title('%d cm' %length[i])

        labels = ['150', '180', '-150', '-120', '-90', '-60']
        xticks = np.array([float(l) for l in labels])
        xticks[xticks < 0] += 360
        labels = [_east_west_formatted(float(l)) for l in labels]

        ax.set_xticklabels(labels)
        ax.set_xticks(xticks)
    plt.suptitle(varname.replace('_', '-'))

    plt.savefig('equatorial_covariance_%s.png' %varname, bbox_inches='tight')
