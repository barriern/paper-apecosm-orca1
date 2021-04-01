import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import sys
sys.path.append('../../nino')
from extract_nino import read_index
import scipy.signal as sig
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import string
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib.ticker as mticker
from cycler import cycler

data = xr.open_dataset('data/ORCA1_JRAC02_CORMSK_CYC1_FINAL_ConstantFields.nc')
data = data.isel(w=[14, 45, 80])
wstep = data['weight_step'].values
print(wstep)

dnino, nino = read_index(filename='data/index/oni.data')
ynino = dnino // 100
iok = np.nonzero((ynino <= 2018) & (ynino >= 1958))
nino = nino[iok]
dnino = dnino[iok]

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree()

#mesh = xr.open_dataset("/Users/Nicolas/Work/sent/apecosm/ORCA1/mesh_mask_eORCA1_v2.2.nc")
mesh = xr.open_dataset("/home/barrier/Work/scientific_communication/articles/ongoing/paper-apecosm-orca1/scripts/data/mesh_mask_eORCA1_v2.2.nc")
lonf = mesh['glamf'].values[0]
latf = mesh['gphif'].values[0]

data = xr.open_dataset("data/eof_oni_annual_density_20.nc")
eof = data['eofmap'].to_masked_array()
pc = data['eofpc'].values
var = data['eofvar'].values * 100
time = data['time'].values
ntime = len(time)

years = np.array([t.year for t in time])
month = np.array([t.month for t in time])
time = np.arange(ntime)
date = ['%.4d-%.2d' %(y, m) for y,m in zip(years, month)]

ncom, nbins, neof = var.shape
coms = ['Epi.']
sizes = ['3cm', '20cm', '90cm']

dictgrid = {'crs':ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
lontext = 120
lattext = 30
            
fig = plt.figure(figsize=(12, 8))
axes_class = (GeoAxes, dict(map_projection=proj))
axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(3, 2), axes_pad=(1.23, 0.4), label_mode='', cbar_mode='each', cbar_pad=0.05)
cbar_axes = axgr.cbar_axes
axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

cpt = 1

letters = list(string.ascii_lowercase)

clim = [50, 5, 2]

for c in range(1):
    for s in range(nbins):
        for e in range(2):

            print("Comm = %d, Size = %d, EOF = %d" %(c, s, e + 1))

            corrcoef = np.corrcoef(nino, pc[c, s, e, :])[0, 1]
            if(corrcoef < 0):
                eof[c, s, e, :, :] *= -1
                pc[c, s, e, : ] *= -1

            eof[c, s, e]  *= wstep[s]

            #ax = plt.subplot(3, 2, cpt, projection=proj)
            ax = axout[cpt - 1]
            temp = np.abs(eof[c, s, e, 1:, 1:])
            perc = np.percentile(temp[temp.mask == False], 98)
            #clim = perc
            
            cs = ax.pcolormesh(lonf, latf, eof[c, s, e, 1:, 1:], transform=proj2, cmap=plt.cm.RdBu_r)
            ax.add_feature(cfeature.LAND, color='lightgray')
            ax.add_feature(cfeature.COASTLINE)
            ax.set_ylim(-40, 40)
            ax.set_xlim(-60, 130)
            cs.set_clim(-clim[s], clim[s])
            title = r'%s, %s, EOF %d (%.2f' %(coms[c], sizes[s], e + 1, var[c, s, e]) + '\%)'
            print(title)
            ax.set_title(title)
            ax.text(lontext, lattext, letters[cpt - 1] + ")", ha='center', va='center', transform=proj, bbox=dicttext)
            #cb = plt.colorbar(cs, orientation='vertical', shrink=1)
            cb = cbar_axes[cpt -1].colorbar(cs)
            #cb.set_label_text('J/m2')
            cb.set_label('J/m2')
            gl = ax.gridlines(**dictgrid)
            gl.xlabels_top = False
            gl.ylabels_right = False
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlocator = mticker.FixedLocator([150, 180, -150, -120, -90, -60])

            cpt += 1
    
            '''
            ax = plt.subplot(3, 2, cpt, xmargin=-0.2)
            print(cpt)
            plt.plot(time, pc[c, s, e, :], color='black')
            plt.fill_between(time, 0, nino, where=nino>0, interpolate=True, color='firebrick')
            plt.fill_between(time, 0, nino, where=nino<0, interpolate=True, color='steelblue')
            ax.set_xlim(time.min(), time.max())
            stride =  5 * 12
            ax.set_xticks(time[::stride])
            ax.set_xticklabels(date[::stride], rotation=45, ha='right')
            ax.set_ylim(-3, 3)
            print(dir(ax))

            cpt += 1

            '''

plt.savefig('eof_full.png', bbox_inches='tight')
plt.close(fig)

'''

            ntime = len(nino) 
            lags = sig.correlation_lags(ntime, ntime)
            ilags = np.nonzero(np.abs(lags) <= 5 * 12)[0]
            corr = sig.correlate(pc[c, s, e, :], nino) / ntime

            lags = lags[ilags]
            corr = corr[ilags]

            imax = np.nonzero(np.abs(corr) == np.abs(corr).max())[0]
            print(lags[imax], corr[imax])
            
            fig = plt.figure()
            plt.plot(lags, corr)
            plt.savefig('correlation_com_%d_size_%s_eof_%d' %(c, s, e), bbox_inches='tight')
            plt.close(fig)
'''
