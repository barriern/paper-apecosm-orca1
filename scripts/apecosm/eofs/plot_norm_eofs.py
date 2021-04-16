import matplotlib
import xarray as xr
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
import matplotlib

signs = -np.ones((100))
index = list(range(84, 93 + 1)) + [99]
index = np.array(index)
signs[index] *= -1

matplotlib.rcParams['image.cmap'] = 'RdBu_r'

const = xr.open_dataset('../../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
length = const['length'].values * 100
wstep = const['weight_step'].values

mesh = xr.open_dataset("../../data/mesh_mask_eORCA1_v2.2.nc")
mesh = mesh.isel(x=slice(58, 228, None), y=slice(140, 233, None))
lonf = mesh['glamf'].values[0]
latf = mesh['gphif'].values[0]

data = xr.open_dataset('data/eof_full_density_20.nc')

eofmap = data['eofmap'].to_masked_array()  #  bins, eof, y, x
eofmap0 = data['eofmap0'].to_masked_array()
eofmap1 = data['eofmap1'].to_masked_array()
eofmap2 = data['eofmap2'].to_masked_array()
eofvar = data['eofvar'].values * 100

eofmap = (eofmap.T * signs).T
eofmap0 = (eofmap0.T * signs).T


nbins, neof, ny, nlat = eofmap.shape

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree()

cc0 = 0.04

with PdfPages('eof_norm.pdf') as pdf:

    for s in range(nbins):

        fig = plt.figure()
        plt.subplots_adjust(hspace=0.3)

        suptitle = 's=%.3d, Length = %.e cm, Variance = %.f' % (s, length[s], eofvar[s, 0]) + '%'
        plt.suptitle(suptitle)

        tp = eofmap[s, 0, 1:, 1:] * wstep[s]
        perc = np.percentile(np.abs(tp[tp.mask == False]), 99.5)

        ax = plt.subplot(2, 1, 1, projection=proj)
        cs = ax.pcolormesh(lonf, latf, tp, transform=proj2)
        cs.set_clim(-perc, perc)
        cb = plt.colorbar(cs, orientation='horizontal')
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
        ax.set_title("Cov.")
        
        tp = eofmap0[s, 0, 1:, 1:]
        perc = np.percentile(np.abs(tp[tp.mask == False]), 99.5)

        ax = plt.subplot(2, 1, 2, projection=proj)
        cs = ax.pcolormesh(lonf, latf, tp, transform=proj2)
        cs.set_clim(-perc, perc)
        cb = plt.colorbar(cs, orientation='horizontal')
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
        ax.set_title("EOF0")
        cs.set_clim(-cc0, cc0)
        
        '''
        tp = eofmap1[s, 0, 1:, 1:]
        perc = np.percentile(np.abs(tp[tp.mask == False]), 99.5)

        ax = plt.subplot(2, 2, 3, projection=proj)
        cs = ax.pcolormesh(lonf, latf, tp, transform=proj2)
        cs.set_clim(-perc, perc)
        cb = plt.colorbar(cs)
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
        ax.set_title("EOF1")
        
        tp = eofmap2[s, 0, 1:, 1:]
        perc = np.percentile(np.abs(tp[tp.mask == False]), 99.5)

        ax = plt.subplot(2, 2, 4, projection=proj)
        cs = ax.pcolormesh(lonf, latf, tp, transform=proj2)
        cs.set_clim(-perc, perc)
        cb = plt.colorbar(cs)
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
        ax.set_title("EOF2")
        '''

        pdf.savefig()
        plt.close(fig)
