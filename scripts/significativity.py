# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import scipy.stats
import numpy as np

student = scipy.stats.t

''' 
Computes the bretherton coefficient
a = autocorrelation at 1-lag for first time-series
b = autocorrelation at 1-lag for second time-series
'''
def bretherthon(a, b):    
    return (1 - a * b) / ( 1 + a * b)

'''
Computes the one-lag autocorrelation
'''
def autocorr(ts1):
    x1 = np.corrcoef(ts1[:-1], ts1[1:])[0, 1]
    return x1

'''
Returns the correlation above which significance is true.
'''
def sig(ts1, ts2, dof=12, proba=0.95):
    x1 = autocorr(ts1)
    x2 = autocorr(ts2)
    print('x1', x1)
    print('x2', x2)
    bres = bretherthon(x1, x2)
    print('bres', bres)
    dl = (len(ts1) - dof) * bres
    print('n = ', len(ts1))
    print('dl', dl)
    tlim = student.interval(proba, dl)[-1]
    print('tlim', tlim)
    output = coef(dl, tlim)
    return output

'''
Computes the correlation based on dl and t value
'''
def coef(dl, tlim):
    return tlim/np.sqrt(dl-2+tlim*tlim)
# -


