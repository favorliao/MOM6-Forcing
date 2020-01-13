#!/bin/bash

# interpolate data into dfs time resolution
cdo inttime,1958-01-01,00:00:00,3hour /tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/origin_slp_era40_1958-01-01to1979-01-01.nc    /tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/regrid_all_slp_era40_1958-01-01to1979-01-01_3hrs.nc
cdo inttime,1979-01-01,00:00:00,3hour /tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/origin_slp_interim_1979-01-01to2018-01-01.nc  /tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/regrid_all_slp_interim_1979-01-01to2018-01-01_3hrs.nc

# interpolate data into dfs time resolution
cdo -f nc -remapcon,/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/data_dfs_mom6/dfs5.2_mom6_u10_1958_10oct2018.nc  /tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/regrid_all_slp_era40_1958-01-01to1979-01-01_3hrs.nc  /tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/regrid_all_slp_era40_1958-01-01to1979-01-01_3hrs_dfsgrid.nc
cdo -f nc -remapcon,/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/data_dfs_mom6/dfs5.2_mom6_u10_1958_10oct2018.nc  /tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/regrid_all_slp_interim_1979-01-01to2018-01-01_3hrs.nc  /tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/regrid_all_slp_interim_1979-01-01to2018-01-01_3hrs_dfsgrid.nc
# split whole files into yearly file
cdo splityear /tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/regrid_all_slp_era40_1958-01-01to1979-01-01_3hrs_dfsgrid.nc     /tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/regrid_yearly_slp_era40_3hrs_dfsgrid
cdo splityear /tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/regrid_all_slp_interim_1979-01-01to2018-01-01_3hrs_dfsgrid.nc   /tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/slp_ecmwf/regrid_yearly_slp_interim_3hrs_dfsgrid

rm regrid_yearly_slp_era40_3hrs_dfsgrid1979.nc # the file is replicated (era40, interim)

# the COBE sst is used instead of ecmwf sst
#cdo inttime,1958-01-01,00:00:00,3hour data_origin/sst_era40_1958-01-01to1979-01-01.nc    data_res/sst_era40_1958-01-01to1979-01-01_3hrs.nc
#cdo inttime,1979-01-01,00:00:00,3hour data_origin/sst_interim_1979-01-01to2018-01-01.nc  data_res/sst_interim_1979-01-01to2018-01-01_3hrs.nc
#cdo -f nc -remapcon,/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/data_dfs/drowned_q2_DFS5.2_y1958.nc  data_res/sst_era40_1958-01-01to1979-01-01_3hrs.nc  ../data_res/sst_era40_1958-01-01to1979-01-01_3hrs_dfsgrid.nc
#cdo -f nc -remapcon,/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/data_dfs/drowned_q2_DFS5.2_y1958.nc  data_res/sst_interim_1979-01-01to2018-01-01_3hrs.nc ../data_res/sst_interim_1979-01-01to2018-01-01_3hrs_dfsgrid.nc
#cdo splityear data_res/sst_era40_1958-01-01to1979-01-01_3hrs_dfsgrid.nc     data_res/sst_era40_3hrs_dfsgrid_yearly
#cdo splityear data_res/sst_interim_1979-01-01to2018-01-01_3hrs_dfsgrid.nc   data_res/sst_interim_3hrs_dfsgrid_yearly
