import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid, ImageGrid

data = xr.open_dataset('full_eof_tpi_oni_corr.nc')
data = data.isel(l=slice(None, -1))
lags = data['lags'].values
length2 = data['l'].values
iok = np.nonzero(length2 <= 100)[0]
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

#corroni = np.abs(corroni)
#corrtpi = np.abs(corrtpi)
cmax = 0.8

print(corrtpi.max())

step = 0.05
lev = np.arange(0, 1 + step, step)
lev = np.arange(-0.2, 0.85 + step, step)

# +
plt.rcParams['font.size'] = 15

fig = plt.figure(figsize=(10, 9), facecolor='white')
axgr = ImageGrid(fig, 111,  nrows_ncols=(2, 1), label_mode='L', aspect=False, 
                 share_all=True, axes_pad=[0.2, 0.4], cbar_mode='single', cbar_size='2%')
cbar_axes = axgr.cbar_axes
axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

ax = axout[0]
cs = ax.contourf(lags, length, corroni, levels=lev)
cl = ax.contour(lags, length, corroni, colors='k', linewidths=0.5, levels=lev)
cl2 = ax.contour(lags, length, corroni, colors='k', linewidths=2, levels=0)
cb = plt.colorbar(cs, cax=cbar_axes[0], shrink=0.1, aspect=0.8)
ax.set_xlabel('Lags (month)')
ax.set_ylabel('Length (cm)')
ax.set_title('ONI')
ax.grid(True)
cb = cbar_axes[0].colorbar(cs)
cb.set_label('Correlation')

ax = axout[1]
cs = ax.contourf(lags, length, corrtpi, levels=lev)
cl = ax.contour(lags, length, corrtpi, colors='k', linewidths=0.5, levels=lev)
cl2 = ax.contour(lags, length, corrtpi, colors='k', linewidths=2, levels=0)
ax.set_xlabel('Lags (month)')
ax.set_ylabel('Length (cm)')
ax.set_title('Filt. TPI')
ax.set_yscale('log')
ax.set_xlim(lags.min(), lags.max())
ax.set_ylim(length.min(), length.max())
ax.grid(True)

plt.savefig('correlations_eof1_oni_tpi.png', bbox_inches='tight')
# -


