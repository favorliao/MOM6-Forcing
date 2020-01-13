#!/bin/bash

# interpolate data into dfs time resolution
cdo mergetime sst.COBESST.*.nc all.sst.COBESST.nc

cdo -f nc -remapcon,/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/data_dfs_mom6/dfs5.2_mom6_u10_1958_10oct2018.nc  all.sst.COBESST.nc  all.sst.COBESST.dfsgrid.nc

cdo inttime,1958-01-01,03:00:00,3hour all.sst.COBESST.dfsgrid.nc   all.sst.COBESST.dfsgrid.3hrs.nc

# interpolate data into dfs time resolution
# split whole files into yearly file
cdo splityear all.sst.COBESST.dfsgrid.3hrs.nc   regrid_3hrs_sst_COBESST_yearly

