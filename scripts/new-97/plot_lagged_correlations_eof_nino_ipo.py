# +
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid, ImageGrid

eof = 0
# -

data = xr.open_dataset('full_eof_tpi_oni_corr_eof_%d.nc' %(eof + 1))
data = data.isel(l=slice(None, -1))
lags = data['lags'].values
length2 = data['l'].values
iok = np.nonzero(length2 <= 100)[0]
nlags = len(lags)

corrtpi = data['corrtpi'].values
corroni = data['corroni'].values
nsize, nlags = corroni.shape
if(eof == 1):
    lagref = 20
else:
    lagref = 0
iii = np.nonzero(lags == lagref)[0][0]
iii

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

fig = plt.figure(figsize=(10, 9), facecolor='white')
axgr = ImageGrid(fig, 111,  nrows_ncols=(2, 1), label_mode='L', aspect=False, 
                 share_all=True, axes_pad=[0.2, 0.4], cbar_mode='single', cbar_size='2%')
cbar_axes = axgr.cbar_axes
axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

ax = axout[0]
#cs = ax.contourf(lags, length, corroni, levels=lev)
cs = ax.pcolor(lags, length, corroni)
cs.set_clim(-cmax, cmax)
cl = ax.contour(lags, length, corroni, colors='k', linewidths=0.5, levels=lev)
cl2 = ax.contour(lags, length, corroni, colors='k', linewidths=2, levels=0)
cb = plt.colorbar(cs, cax=cbar_axes[0], shrink=0.1, aspect=0.8)
ax.set_xlabel('Lags (month)')
ax.set_ylabel('Length (cm)')
ax.set_title('ONI')
ax.grid(True)
cb = cbar_axes[0].colorbar(cs)
cb.set_label('Correlation PC %d' %(eof + 1))
ax.axvline(lagref)

ax = axout[1]
cs = ax.contourf(lags, length, corrtpi, levels=lev)
cs = ax.pcolor(lags, length, corrtpi)
cs.set_clim(-cmax, cmax)
cl = ax.contour(lags, length, corrtpi, colors='k', linewidths=0.5, levels=lev)
cl2 = ax.contour(lags, length, corrtpi, colors='k', linewidths=2, levels=0)
ax.set_xlabel('Lags (month)')
ax.set_ylabel('Length (cm)')
ax.set_title('Filt. TPI')
ax.set_yscale('log')
ax.set_xlim(lags.min(), lags.max())
ax.set_ylim(length.min(), length.max())
ax.grid(True)

plt.savefig('correlations_eof_oni_tpi_eof_%d.png' %(eof + 1), bbox_inches='tight')
# -
plt.plot(signs_oni)




