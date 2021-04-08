import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import scipy.signal as sig
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import string
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib.ticker as mticker
from cycler import cycler
import apecosm.ts as ts

dirin = '/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/JRA_CO2/0-200'

ilat = slice(128, 246, None)
ilon = slice(58, 228, None)

mesh = xr.open_dataset("/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/corrected_mesh_mask_eORCA1_v2.2.nc")
mesh = mesh.isel(y=ilat, x=ilon)
lonf = mesh['glamf'].values[0]
latf = mesh['gphif'].values[0]

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree()
    
dataset = xr.open_mfdataset('%s/*.nc' %(dirin))
dataset = dataset.isel(y=ilat, x=ilon)  # (ntime, y, x, w=3)

def plot_seasonal_cycle(varname):

    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@", varname)

    data = dataset[varname].to_masked_array()
    print(data.shape)
    mask = np.ma.getmaskarray(data)[:12, :, :]
    data, anom = ts.get_monthly_clim(data)
    print(data.shape)
    del(anom)
    data = np.ma.masked_where(mask == True, data)
    ntime, ny, nx = data.shape

    dictgrid = {'crs':ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}
    dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)

    fig = plt.figure(figsize=(12, 8))
    axes_class = (GeoAxes, dict(map_projection=proj))
    axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(4, 3), cbar_mode='single', cbar_pad=0.05, label_mode='', axes_pad=(0.1, 0.4))
    cbar_axes = axgr.cbar_axes
    axout = list(enumerate(axgr))
    axout = [p[1] for p in axout]

    temp = data
    iok = np.nonzero(np.ma.getmaskarray(temp) == False)
    off = 1
    cmin = np.percentile(temp[iok], off)
    cmax = np.percentile(temp[iok], 100 - off)

    if (('uo' in varname) | ('vo' in varname)):
        cm = max(-cmin, cmax)
        cmin = -cm
        cmax = cm

    cpt = 1
    for t in range(12):

        ax = axout[cpt - 1]

        print(temp[t, 1:, 1:].min(), temp[t, 1:, 1:].max())
        
        cs = ax.pcolormesh(lonf, latf, temp[t, 1:, 1:], transform=proj2, cmap=plt.cm.jet)
        ax.add_feature(cfeature.LAND, color='lightgray')
        ax.add_feature(cfeature.COASTLINE)
        ax.set_ylim(-40, 40)
        ax.set_xlim(-60, 130)
        cs.set_clim(cmin, cmax)
        #ax.set_title('Month=%.2d, %s, %s, %s' %(t + 1, coms[0], sizes[s], varname.replace('_', '-'))) 
        ax.set_title('Month = %.2d' %(t + 1))
        #cb = cbar_axes[cpt - 1].colorbar(cs)
        cb = cbar_axes[0].colorbar(cs)
        
        cpt += 1

    plt.savefig('clim_%s.png' %(varname), bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':

    varlist = ['thetao', 'uo', 'vo', 'ZOO', 'ZOO2', 'PHY2', 'O2', 'GOC']
    for v in varlist:
        plot_seasonal_cycle(v)
