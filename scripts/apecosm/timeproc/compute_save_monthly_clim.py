import xarray as xr
import numpy as np

dirin = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output/'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data/'

def compute_clim(varname):

    print('++++++++++++++++++++ Processing ', varname)

    data = xr.open_mfdataset("%s/*%s*nc" %(dirin, varname))
    print(data)
    ntime = data.dims['time']
    nyears = ntime // 12

    dims = data[varname].dims

    index = np.arange(12).astype(np.int)

    for y in range(nyears):
        
        if(y == 0):
            output = data[varname].isel(time=index).values
        else:
            output += data[varname].isel(time=index).values

        print(index)

        index += 12

    output /= nyears

    dsout = xr.Dataset()
    dsout[varname] = (dims, output)
    dsout.to_netcdf("%s/clim_%s.nc" %(dirout, varname))

if __name__ == '__main__':

    varlist = ['u_active', 'v_active', 'u_passive', 'v_passive', 'zdiff_trend', 'mdiff_trend', 'madv_trend', 'zadv_trend']
    for v in varlist:
        compute_clim(v)
