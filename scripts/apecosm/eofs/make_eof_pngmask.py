import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

data = xr.open_dataset("../../data/mesh_mask_eORCA1_v2.2.nc")
data = data.isel(t=0, z=0)

tmask = data['tmask'].values
lon = data['glamt'].values
lat = data['gphit'].values
lonf = data['glamf'].values
latf = data['gphif'].values

latmax = 20

test = (np.abs(lat) <= latmax)
test2 = (lon >= 130) | (lon <= -60)
test = test & test2

tmask[test] = tmask[test] * 3

plt.figure()
plt.imsave('temp_mask.png', tmask)




