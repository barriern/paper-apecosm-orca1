# +
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid, ImageGrid

eof = 0
latmax = 20
lonmin = 150
lonmax = -120
# -

data = xr.open_dataset('full_eof_tpi_oni_corr_eof_%d.nc' %(eof + 1))
data = data.isel(l=slice(None, -1))
lags = data['lags'].values
length2 = data['l'].values
iok = np.nonzero(length2 <= 100)[0]
nlags = len(lags)

# +
filename = 'data/full_eof_pacific_OOPE_latmax_%d_lonmin_%d_lonmax_%d.nc' %(latmax, lonmin, lonmax)
filename

eofvar = xr.open_dataset(filename)
eofvar = eofvar['eofvar'].isel(w=slice(None, -1))
eofvar = eofvar.rename({'w': 'l'})
eofvar['l'] = length2
eofvar = eofvar[:, eof]
eofvar
# -

corrtpi = data['corrtpi'].values  # nsize, nlags
corroni = data['corroni'].values
nsize, nlags = corroni.shape
if(eof == 1):
    lagref = 20
else:
    lagref = 0
iii = np.nonzero(lags == lagref)[0][0]
iii
corroni.shape

signs_oni = np.sign(corroni[:, iii])  # nsize
signs_oni.shape
corroni[signs_oni < 0, :] *= -1
corrtpi[signs_oni < 0, :] *= -1

il = np.nonzero(length2 <= 100)[0]
length = length2[il]
corrtpi = corrtpi[il, :]
corroni = corroni[il, :]

ilags = np.nonzero(lags >= 0)[0]
corrtpi = corrtpi[:, ilags]
corroni = corroni[:, ilags]
lags = lags[ilags]

cmax = 0.8

step = 0.05
lev = np.arange(0, 1 + step, step)
lev = np.arange(-0.2, 0.85 + step, step)
lev = np.arange(-cmax, cmax + step, step)

# +
plt.rcParams['font.size'] = 15

fig = plt.figure(figsize=(10, 14), facecolor='white')
#axgr = ImageGrid(fig, 311,  nrows_ncols=(2, 1), label_mode='L', aspect=False, 
#                 share_all=False, axes_pad=[0.2, 0.4], cbar_mode='single', cbar_size='5%', cbar_location='bottom', cbar_pad='0.1%')
#cbar_axes = axgr.cbar_axes
#axout = list(enumerate(axgr))
#axout = [p[1] for p in axout]

ax = plt.subplot(311)
#cs = ax.contourf(lags, length, corroni, levels=lev)
cs = ax.pcolor(length, lags, corroni.T, shading='auto')
cs.set_clim(-cmax, cmax)
cl = ax.contour(length, lags,  corroni.T, colors='k', linewidths=0.5, levels=lev)
cl2 = ax.contour(length, lags, corroni.T, colors='k', linewidths=2, levels=0)
cb = plt.colorbar(cs)
#ax.set_xlabel('Lags (month)')
ax.set_ylabel('Length (cm)')
ax.set_title('ONI')
ax.grid(True)
#cb = cbar_axes[0].colorbar(cs)
cb.set_label('Correlation PC %d' %(eof + 1))
#ax.axvline(lagref)
plt.tick_params(labelbottom=False)

ax = plt.subplot(312, sharex=ax)
#cs = ax.contourf(lags, length, corrtpi, levels=lev)
cs = ax.pcolor(length, lags, corrtpi.T, shading='auto')
cs.set_clim(-cmax, cmax)
cl = ax.contour(length, lags,  corrtpi.T, colors='k', linewidths=0.5, levels=lev)
cl2 = ax.contour(length, lags, corrtpi.T, colors='k', linewidths=2, levels=0)
#cb = plt.colorbar(cs, cax=cbar_axes[1], shrink=0.1, aspect=0.8)
cb = plt.colorbar(cs)
ax.set_ylabel('Lags (month)')
ax.set_title('Filt. TPI')
cb.set_label('Correlation PC %d' %(eof + 1))
ax.set_ylim(lags.min(), lags.max())
ax.set_xlim(length.min(), length.max())
ax.grid(True)
ax.set_xscale('log')
plt.tick_params(labelbottom=False)

pos1 = ax.get_position() # get the original position 
offset = 0.25
pos2 = [pos1.x0, pos1.y0 -offset,  pos1.width, pos1.height]

ax = plt.axes(pos2)
ax.plot(length, eofvar[il].values)
ax.set_xlabel('Length (cm)')
ax.set_xscale('log')
ax.set_xlim(length.min(), length.max())
ax.set_ylabel('\%')
ax.set_title('Explained variance')
plt.grid(True)
plt.savefig('correlations_eof_oni_tpi_eof_%d.png' %(eof + 1))