import numpy as np

def find_percentile(data, percentage=1):

    ''' 
    Extract percentile to saturate the colormaps.
    They are computed from unmasked arrays

    :param numpy.array data: Data array
    :param float percentage: Percentage used to
     saturate the colormap.

    :return: A tuple containing the lower and upper bounds (cmin, cmax)

    '''
    data = np.ma.masked_where(np.isnan(data), data)
    iok = np.nonzero(np.logical_not(np.ma.getmaskarray(data)))
    temp = data[iok]

    cmin = np.percentile(np.ravel(temp), percentage)
    cmax = np.percentile(np.ravel(temp), 100 - percentage)

    return cmin, cmax

def get_monthly_clim(var):
    
    ntime = var.shape[0]

    nyears = ntime // 12

    index = np.arange(12)

    for i in range(nyears):
        if(i == 0):
            clim = var[index]
        else:
            clim += var[index]
        index += 12
    clim /= nyears

    anom = np.zeros(var.shape)
    index = np.arange(12)
    for i in range(nyears):
        anom[index] = var[index] - clim
        index += 12

    return clim, anom
