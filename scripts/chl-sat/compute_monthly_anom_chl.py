from glob import glob
import numpy as np
import xarray as xr

dataclim = xr.open_dataset("ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-CLIM-fv5.0.nc")
clim = dataclim['clim_chlor_a'].to_masked_array()

filelist = np.sort(glob('ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-[0-9]*-fv5.0.nc'))

for f in filelist:

    data = xr.open_dataset(f)
    month = data['time.month'].values[0] - 1 
    print(f, month)

    chl = data['chlor_a']

    chl = chl - clim[month:month + 1]

    fout = 'anom_%s' %f

    chl.to_netcdf(fout)
