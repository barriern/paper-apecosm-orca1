import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid, AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

constant = xr.open_dataset('/home/barrier/Work/apecosm/ORCA1/figures_orca1/plot_hov/ORCA1_JRA_CO2_CYC4_ConstantFields.nc')
wstep = constant['weight_step'].values

index = slice(None, None)
data = xr.open_dataset('/home/barrier/Work/scientific_communication/articles/ongoing/paper-apecosm-orca1/scripts/data/ORCA1_JRA_CO2_CYC4_ConstantFields.nc')
length = data['length'].values[index] * 100
wstep = data['weight_step'].values[index]

data = xr.open_dataset('/home/barrier/Work/apecosm/ORCA1/DATA_APE_ORCA1/apecosm/data/corrected_mesh_mask_eORCA1_v2.2.nc')
data = data.isel(t=0)

lonf = data['glamf'].values
latf = data['gphif'].values
lon = data['glamt'].values
lat = data['gphit'].values
tmask = data['tmask'].values[0]
nlat, nlon = tmask.shape

lonbis = lon.copy()
lonbis[lonbis < 0] += 360

test = (lonbis >= 130) & (lonbis <= 300)
test = test & (np.abs(lat) <=40)
test = test & (tmask > 0)

ilat, ilon = np.nonzero(test == True)

data = xr.open_dataset('lagged_covariance_monthly_epi_allsizes.nc')
data = data.sel(lags=slice(0, 60))
data = data.isel(w=index)
cov = data['cov'].to_masked_array()
lags = data['lags'].values

cov = cov * wstep[:, np.newaxis, np.newaxis]
cov /= (61 * 12 - 1)

nw, npoints, nlags = cov.shape

perc = np.percentile(np.abs(cov), 99.9, axis=(1, 2))

projection2 = ccrs.PlateCarree()
projection = ccrs.PlateCarree(central_longitude=180)
axes_class = (GeoAxes, dict(map_projection=projection))

for s in range(nw):

    covtemp = cov[s]
    ltemp = length[s]

    for p in range(nlags):

        print('Processing lag %d / %d' %(p, nlags))

        fig = plt.figure()
        ax = plt.axes(projection=projection)
        temp = np.zeros((nlat, nlon))
        temp[ilat, ilon] = covtemp[:, p]
        temp = np.ma.masked_where(temp == 0, temp)
        cs = ax.pcolormesh(lonf, latf, temp[1:, 1:], transform=projection2, cmap=plt.cm.get_cmap('RdBu_r'))
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.LAND)
        cb = plt.colorbar(cs, orientation='horizontal')
        cb.set_label('Covariance OOPE/ONI')
        plt.title('Length = %.2e, Lag = %.2d months' %(ltemp, p))
        cs.set_clim(-perc[s], perc[s])
        ax.set_ylim(-40, 40)
        ax.set_xlim(-60, 130)
        plt.savefig('figs/length_%.3d_lag_%.3d' %(s, p + 1), bbox_inches='tight')
        plt.close(fig)
