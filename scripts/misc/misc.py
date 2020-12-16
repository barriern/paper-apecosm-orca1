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
