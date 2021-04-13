import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid, ImageGrid

data = xr.open_dataset('full_eof_tpi_oni_corr.nc')
data = data.isel(l=slice(None, -1))
lags = data['lags'].values
length2 = data['l'].values
nlags = len(lags)

corrtpi = data['corrtpi'].values  # size, lags
corroni = data['corroni'].values  # size, lags

il = np.nonzero(length2 <= 100)[0]
length = length2[il]
corrtpi = corrtpi[il, :]
corroni = corroni[il, :]

ilags = np.nonzero(lags >= 0)[0]
corrtpi = corrtpi[:, ilags]
corroni = corroni[:, ilags]
lags = lags[ilags]

isize = np.array([ 14, 45, 80])

corroni = np.abs(corroni)
corrtpi = np.abs(corrtpi)
cmax = 0.8

print(corrtpi.max())

step = 0.05
lev = np.arange(0, 1 + step, step)

fig = plt.figure()
axgr = ImageGrid(fig, 111,  nrows_ncols=(2, 1), label_mode='L', aspect=False, share_all=True, axes_pad=[0.2, 0.3], cbar_mode='single')
cbar_axes = axgr.cbar_axes
axout = list(enumerate(axgr))
axout = [p[1] for p in axout]
print(axout)
print(cbar_axes)

ax = axout[0]
cs = ax.pcolormesh(lags, length, corroni, cmap=plt.cm.jet)
cl = ax.contour(lags, length, corroni, colors='k', linewidths=0.5, levels=lev)
cbar_axes[0].colorbar(cs)
cs.set_clim(0, cmax)
ax.set_xlabel('Lags (month)')
ax.set_ylabel('Length (cm)')
ax.set_title('ONI')
#for l in length2[isize]:
#    ax.axhline(l, linestyle='--', color='k')

ax = axout[1]
cs = ax.pcolormesh(lags, length, corrtpi, cmap=plt.cm.jet)
cl = ax.contour(lags, length, corrtpi, colors='k', linewidths=0.5, levels=lev)
cb = cbar_axes[1].colorbar(cs)
cs.set_clim(0, cmax)
ax.set_xlabel('Lags (month)')
ax.set_ylabel('Length (cm)')
ax.set_title('Filt. TPI')
#for l in length2[isize]:
#    ax.axhline(l, linestyle='--', color='k')
ax.set_yscale('log')
ax.set_xlim(lags.min(), lags.max())
ax.set_ylim(length.min(), length.max())
try:
    cb.set_label('Correlation')
except:
    cb.set_label_text('Correlation')

plt.savefig('correlations_eof1_oni_tpi.png', bbox_inches='tight')



