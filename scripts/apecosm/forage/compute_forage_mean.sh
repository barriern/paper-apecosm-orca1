#!/bin/bash

prefix=yearly_mean_forage
prefix=meridional_-150_forage

cd /home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data
ncra -O  ${prefix}_year_*.nc ${prefix}.nc
