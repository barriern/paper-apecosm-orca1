from glob import glob
import xarray as xr
import numpy as np

dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output/'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data/'
dirout = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output/tsanalysis'

ilat = slice(128, 246, None)
ilon = slice(58, 228, None)

def prepare_for_timeseries(varname):

    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ processing ", varname)

    filelist = np.sort(glob('%s/*%s*nc' %(dirin, varname)))
    for f in filelist[:]:

        print("--- ", f)

        data = xr.open_dataset(f)
        data = data.isel(y=ilat, x=ilon, community=0)
        data = data[varname].T  # w, lon, lat, time
        dims = data.dims
        data = data.values
        if(f == filelist[0]):
            output = data
        else:
            output = np.concatenate((output, data), axis=-1)
    dsout = xr.Dataset()
    dsout[varname] = (dims, output)
    dsout.to_netcdf('%s/flipped_%s' %(dirout, varname))

if __name__ == '__main__':

    prepare_for_timeseries('OOPE')
    prepare_for_timeseries('repfonct_day')
