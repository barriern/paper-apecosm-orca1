import matplotlib
import xarray as xr
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

const = xr.open_dataset('../../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
length = const['length'].values * 100
wstep = const['weight_step'].values

mesh = xr.open_dataset("../../data/mesh_mask_eORCA1_v2.2.nc")
mesh = mesh.isel(x=slice(58, 228, None), y=slice(140, 233, None))
lonf = mesh['glamf'].values[0]
latf = mesh['gphif'].values[0]

data = xr.open_dataset('data/eof_full_density_20.nc')

eofmap = data['eofmap'].to_masked_array()  #  bins, eof, y, x
eofmap2 = data['eofmap2'].to_masked_array()
eofvar = data['eofvar'].values * 100


nbins, neof, ny, nlat = eofmap.shape

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree()

with PdfPages('eof_norm.pdf') as pdf:

    for s in range(nbins):

        fig = plt.figure()

        suptitle = 'Length = %.e cm, Variance = %.f' % (length[s], eofvar[s, 0])
        plt.suptitle(suptitle)

        tp = eofmap[s, 0, 1:, 1:] * wstep[s]
        perc = np.percentile(np.abs(tp[tp.mask == False]), 99.5)

        ax = plt.subplot(2, 1, 1, projection=proj)
        cs = ax.pcolormesh(lonf, latf, tp, transform=proj2)
        cs.set_clim(-perc, perc)
        cb = plt.colorbar(cs)
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
        
        tp = eofmap2[s, 0, 1:, 1:]
        perc = np.percentile(np.abs(tp[tp.mask == False]), 99.5)

        ax = plt.subplot(2, 1, 2, projection=proj)
        cs = ax.pcolormesh(lonf, latf, tp, transform=proj2)
        cs.set_clim(-perc, perc)
        cb = plt.colorbar(cs)
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)

        pdf.savefig()

        break
        
