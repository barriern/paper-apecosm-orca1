import xarray as xr
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("../../nino/")
from extract_nino import read_index
from cycler import cycler
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['axes.prop_cycle'] = cycler('color', ['black', 'green', 'orange', 'brown'])

dnino, nino = read_index(filename='../../data/index/oni.data')
ynino = dnino // 100
iok = np.nonzero((ynino <= 2018) & (ynino >= 1958))
nino = nino[iok]
dnino = dnino[iok]

ieof = 1

ilon = slice(58, 229, None)
ilat = slice(114, 265, None)

const = xr.open_dataset("../../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc")
length = const['length'].values * 100
weight_step = const['weight_step'].values

mesh = xr.open_dataset("../../data/mesh_mask_eORCA1_v2.2.nc")
mesh = mesh.isel(x=ilon, y=ilat)
lonf = mesh['glamf'].values[0]
latf = mesh['gphif'].values[0]

data2 = xr.open_dataset("data/eof_full_density_20.nc")
data2 = data2.isel(x=slice(1, None, None), y=slice(1, None, None), eof=ieof)

data3 = xr.open_dataset("data/eof_full_density_30.nc")
data3 = data3.isel(x=slice(1, None, None), y=slice(1, None, None), eof=ieof)

data4 = xr.open_dataset("data/eof_full_density_40.nc")
data4 = data4.isel(x=slice(1, None, None), y=slice(1, None, None), eof=ieof)

data5 = xr.open_dataset("data/eof_full_density_50.nc")
data5 = data5.isel(x=slice(1, None, None), y=slice(1, None, None), eof=ieof)

time = data5['time'].values
year = [t.year for t in time]
month = [t.month for t in time]
labels = ['%.4d-%.2d' %(y, m) for y, m in zip(year, month)]
time = np.arange(0, len(time))

# eofmap: bins, eof, y,x 
# eofpc: bins, eof, time

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree()
shrink = 0.8

clim = [50, 5, 2]
if(ieof == 1):
    clim = clim

with PdfPages('eof%d_all_lat.pdf' %(ieof + 1)) as pdf:

    toto = 0

    for s in [14, 45, 80]:

        ccc = clim[toto]

        pc2 = data2['eofpc'].isel(bins=s).to_masked_array()
        pc3 = data3['eofpc'].isel(bins=s).to_masked_array()
        pc4 = data4['eofpc'].isel(bins=s).to_masked_array()
        pc5 = data5['eofpc'].isel(bins=s).to_masked_array()

        # calculation of correlations
        corr2 = np.corrcoef(pc2, nino)[0, 1]
        corr3 = np.corrcoef(pc3, nino)[0, 1]
        corr4 = np.corrcoef(pc4, nino)[0, 1]
        corr5 = np.corrcoef(pc5, nino)[0, 1]

        fig = plt.figure()
        sign = np.sign(corr2)
        corr2 = sign * corr2
        plt.suptitle('Length = %.f cm' %(length[s]))
        cpt = 1
        ax = plt.subplot(2, 2, cpt, projection=proj)
        cs = ax.pcolormesh(lonf, latf, sign * data2['eofmap'].isel(bins=s).to_masked_array() * weight_step[s], transform=proj2)
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
        plt.colorbar(cs, orientation='horizontal', shrink=shrink)
        cs.set_clim(-ccc, ccc)
        ax.set_title('Corr = %.2f' %corr2)
        
        sign = np.sign(corr3)
        corr3 = sign * corr3
        cpt += 1
        ax = plt.subplot(2, 2, cpt, projection=proj)
        cs = ax.pcolormesh(lonf, latf, sign * data3['eofmap'].isel(bins=s).to_masked_array() * weight_step[s], transform=proj2)
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
        plt.colorbar(cs, orientation='horizontal', shrink=shrink)
        cs.set_clim(-ccc, ccc)
        ax.set_title('Corr = %.2f' %corr3)
        
        cpt += 1
        ax = plt.subplot(2, 2, cpt, projection=proj)
        sign = np.sign(corr4)
        corr4 = sign * corr4
        cs = ax.pcolormesh(lonf, latf, sign * data4['eofmap'].isel(bins=s).to_masked_array() * weight_step[s], transform=proj2)
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
        plt.colorbar(cs, orientation='horizontal', shrink=shrink)
        cs.set_clim(-ccc, ccc)
        ax.set_title('Corr = %.2f' %corr4)
        
        cpt += 1
        ax = plt.subplot(2, 2, cpt, projection=proj)
        sign = np.sign(corr5)
        corr5 = sign * corr5
        cs = ax.pcolormesh(lonf, latf, sign * data5['eofmap'].isel(bins=s).to_masked_array() * weight_step[s], transform=proj2)
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
        plt.colorbar(cs, orientation='horizontal', shrink=shrink)
        cs.set_clim(-ccc, ccc)
        ax.set_title('Corr = %.2f' %corr5)

        toto += 1
        
        pdf.savefig()
        plt.close(fig)

with PdfPages('pc%d_all_lat.pdf' %(ieof + 1)) as pdf:

    for s in [14, 45, 80]:
        fig = plt.figure()
        
        plt.title('Length = %.f cm' %(length[s]) )
        
        pc2 = data2['eofpc'].isel(bins=s).to_masked_array()
        pc3 = data3['eofpc'].isel(bins=s).to_masked_array()
        pc4 = data4['eofpc'].isel(bins=s).to_masked_array()
        pc5 = data5['eofpc'].isel(bins=s).to_masked_array()

        # calculation of correlations
        corr2 = np.corrcoef(pc2, nino)[0, 1]
        corr3 = np.corrcoef(pc3, nino)[0, 1]
        corr4 = np.corrcoef(pc4, nino)[0, 1]
        corr5 = np.corrcoef(pc5, nino)[0, 1]

        sign2 = np.sign(corr2)
        sign3 = np.sign(corr3)
        sign4 = np.sign(corr4)
        sign5 = np.sign(corr5)

        ax = plt.gca()
        alpha = 0.5
        ax.fill_between(time, 0, nino, where=(nino>0), interpolate=True, color='firebrick', alpha=alpha)
        ax.fill_between(time, 0, nino, where=(nino<0), interpolate=True, color='steelblue', alpha=alpha)
        ax.plot(time, sign2 * data2['eofpc'].isel(bins=s).to_masked_array(), label='20')
        ax.plot(time, sign3 * data3['eofpc'].isel(bins=s).to_masked_array(), label='30')
        ax.plot(time, sign4 * data4['eofpc'].isel(bins=s).to_masked_array(), label='40')
        ax.plot(time, sign5 * data5['eofpc'].isel(bins=s).to_masked_array(), label='50')
        ax.set_ylim(-3, 3)
        plt.grid()
        stride = 3 * 12
        plt.legend(ncol=2, fontsize=8)
        ax.set_xticks(time[::stride])
        ax.set_xticklabels(labels[::stride], ha='right', rotation=45)
        
        pdf.savefig()
        plt.close(fig)
