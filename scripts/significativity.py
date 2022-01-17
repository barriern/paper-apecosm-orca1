# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import scipy.stats

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
    x1 = np.corrcoef(ts1[:-1], ts1[1:])
    return x1

'''
Returns the correlation above which significance is true.
'''
def sig(ts1, ts2, dof=12):
    x1 = autocorr(ts1)
    x2 = autocorr(ts2)
    bres = bretherton(x1, x2)
    dl = (len(ts1) - dof) * bres
    tlim = student.interval(proba, dl)
    output = coef(dl, tlim)

'''
Computes the correlation based on dl and t value
'''
def coef(dl, tlim):
    return tlim/np.sqrt(dl-2+tlim*tlim)
