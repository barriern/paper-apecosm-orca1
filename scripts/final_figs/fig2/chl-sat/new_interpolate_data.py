import xarray as xr
from glob import glob
import numpy as np
import os.path

filelist = np.sort(glob('raw/ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-[0-9]*nc'))

latout = np.arange(-89.5, 89.5 + 1, 1)[::-1]
lonout = np.arange(0.5, 360, 1)

nlatout = len(latout)
nlonout = len(lonout)

dataout = xr.Dataset({'lat': latout, 'lon':lonout})

for f in filelist:

    print(f)

    datain = xr.open_dataset(f)

    # conversion from Atl cent. to Pac. cent.
    lonin = datain['lon']
    lonin = (lonin + 360) % 360
    datain['lon'] = lonin
    datain = datain.sortby(datain['lon'])

    time = datain['time']
    chlin = datain['chlor_a'].to_masked_array()

    latin = datain['lat'].values

    # init. output and weights arrays
    chlout = np.zeros((1, nlatout, nlonout))
    countout = np.zeros((1, nlatout, nlonout), dtype=np.int)
    weights = np.cos(np.deg2rad(latin))
    
    # slice for latitude
    ilat = np.array([0, 24], dtype=np.int)
    for j in range(nlatout):
        
        # converts  sl array to sl object
        sllat = slice(ilat[0], ilat[1])

        # slice for longitude
        ilon = np.array([0, 24], dtype=np.int)

        for i in range(nlonout):
            # converts  sl array to sl object
            sllon  = slice(ilon[0], ilon[1])
            
            chltemp = chlin[:, sllat, sllon] # 1(time), 24(lat), 24(lon)
            wtemp = weights[sllat]  # 24(lat)
            iii, jjj, kkk = np.nonzero(np.ma.getmaskarray(chltemp) == False)
            npoints = len(iii)   # number of points non masked

            if(npoints != 0):
                num = chltemp[iii, jjj, kkk] * wtemp[jjj]  # npoints
                den = wtemp[jjj]  # npoints
                countout[:, j, i] = npoints * 100 / np.prod(chltemp.shape)
                chlout[:, j, i] = np.sum(num) / np.sum(den)
            
            ilon += 24
        ilat += 24
    
    fout = 'interpolated_%s' %(os.path.basename(f))
    dataout = xr.Dataset({'lon': lonout, 'lat': latout})
    dataout['time'] = time
    dataout['chlor_a'] = (['time', 'lat', 'lon'], chlout)
    dataout['count_a'] = (['time', 'lat', 'lon'], countout)
    dataout.to_netcdf(fout, unlimited_dims='time')
