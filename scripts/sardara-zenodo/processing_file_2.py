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

# # Processing Sardara #2
#
# The aim here is to convert date strings into dates and to remove data which are not monthly. Then, a `year` and a `month` columns is created.

# +
import pandas as pd
import numpy as np
from datetime import datetime

data = pd.read_csv('processed_sardara_1.csv', index_col=0)
data
# -

datestr_start = data['time_start'].values
datestr_start

datestr_end = data['time_end'].values
datestr_end

from datetime import datetime
def process_dates(datestr):
    output = []
    for d in datestr:
        year, month, day = np.array(d[:10].split('-')).astype(int)
        output.append(datetime(year, month, day))
    return np.array(output)


date_end = process_dates(datestr_end)
date_start =  process_dates(datestr_start)

delta = np.array([(end - start).days for end, start in zip(date_end, date_start)])
bool_delta = (delta <= 31)

iok = np.nonzero(bool_delta)[0]
dataout = data.iloc[iok, :]
dataout

years = [d.year for d in date_end[iok]]
months = [d.month for d in date_end[iok]]

dataout = dataout.assign(year = years)

dataout = dataout.assign(month = months)

dataout.to_csv('processed_sardara_2.csv')


