#!/bin/bash

ori_dir='/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/data_dfs_mom6'
mer_dir='/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/res_dfs_mom6_all'

#declare -a varall=("precip" "q10" "radlw" "radsw" "slp" "snow" "t10" "u10" "v10")
#declare -a varall=("precip" "radlw" "radsw" "slp" "u10" "v10" "q10" "t10")
declare -a varall=("snow")

varlength=${#varall[@]}
# here to select the variable for merging
for (( i=1; i<${varlength}+1; i++ ));
#for (( i=1; i<3; i++ ));
do
var=${varall[$i-1]}
file_out="$mer_dir/$var.DFS.all.1958-2015.nc"
file_out_float="$mer_dir/$var.DFS.all.1958-2015.float2.nc"

   echo "$i/$varlength merge $var ..."
   echo "output file name: $file_out_float"
   cdo -b f32 mergetime $ori_dir/*$var*.nc $file_out_float
   echo "     change time att $var ..."
   ncatted -O -a modulo,TIME,c,c,' ' -a modulo_beg,TIME,c,c,'1958-01-01 00:00:00' -a modulo_end,TIME,c,c,'2016-01-01 00:00:00' $file_out_float
#   cdo -b f32 -f nc copy $file_out $file_out_float
   echo "     set calendar $var ..."
   ncatted  -a  calendar,TIME,o,c,"gregorian" $file_out_float
done   


