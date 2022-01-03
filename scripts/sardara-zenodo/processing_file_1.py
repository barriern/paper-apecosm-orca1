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

# # Processing Sardara #1
#
# The present script aims at reducing the size of the Sardara file. The first processing is to reduce the size of the input file by extracting the SKJ and YFT species captured by purse seine (PS) gears.

import pandas as pd
import numpy as np

data = pd.read_csv('global_catch_tunaatlasird_level2.csv')
data

species = data['species'].values
species

gear = data['gear_group'].values
gear

bool_species = (species == 'SKJ') | (species == 'YFT')
bool_species

bool_gear = (gear == 'PS')

iok = np.nonzero(bool_gear & bool_species)[0]
iok

dataout = data.iloc[iok, :]
dataout

dataout.to_csv('processed_sardara_1.csv')
