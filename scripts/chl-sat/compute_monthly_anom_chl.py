import os.path
from glob import glob
import numpy as np
import xarray as xr

dirin = '/home1/datawork/nbarrier/chl-data/'

dataclim = xr.open_dataset("%s/ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-CLIM-fv5.0.nc" %dirin)
clim = dataclim['clim_chlor_a'].to_masked_array()

filelist = np.sort(glob('%s/ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-[0-9]*-fv5.0.nc' %dirin))

for f in filelist:

    data = xr.open_dataset(f)
    month = data['time.month'].values[0] - 1 

    chl = data['chlor_a']

    chl = chl - clim[month:month + 1]

    fout = os.path.basename(f)
    fout = '%s/anom_%s' %(dirin, fout)

    chl.to_netcdf(fout)
