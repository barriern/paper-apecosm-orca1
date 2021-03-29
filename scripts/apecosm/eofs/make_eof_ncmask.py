import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os.path

data = xr.open_dataset("../../data/mesh_mask_eORCA1_v2.2.nc")
data = data.isel(t=0, z=0)

tmask = data['tmask'].values
lon = data['glamt'].values
lat = data['gphit'].values
lonf = data['glamf'].values
latf = data['gphif'].values

data = plt.imread('temp_mask.png')[:, :, 0]
data = data[::-1, :]

unique = np.unique(data)
index = np.nonzero(data == unique[1])

tmask = np.zeros(data.shape)
tmask[index] = 1

data = xr.Dataset()
data['mask'] = (['y', 'x'], tmask)
data.attrs['file'] = os.path.realpath(__file__)
data.to_netcdf('eof_mask.nc')

