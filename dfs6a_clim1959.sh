#!/bin/bash

module load cdo/1.8.2

ori_dir='/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/res_dfs_mom6_3yrs'
res_dir='/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/res_dfs_mom6_clim/res_step1sh'
#res_dir='/tigress/GEOCLIM/LRGROUP/datasets/JRA55-do-v1.3merge/res_clim_detrend'

declare -a varall=("precip" "q10" "radlw" "radsw" "slp" "snow" "t10" "u10" "v10")

cd $res_dir
varlength=${#varall[@]}
# here to select the variable for merging
for (( i=1; i<${varlength}+1; i++ ));
#for (( i=5; i<=5; i++ ));
do
   var=${varall[$i-1]}
   file_out="$res_dir/$var.DFS.clim.1959.nc"
   file_in=(`ls $ori_dir/*$var*1959*.nc`) 
   echo "$i/$varlength clim 1959 $var ..."
   echo "input file name: $file_in"
   echo "output file name1: $file_out"

   if [[ ! -e "$file_out" ]]; then
      cdo -O -selyear,1959 $file_in $file_out
   fi
   echo "     change time att ..."
   ncatted -a modulo,TIME,c,c,' ' -a modulo_beg,TIME,c,c,"1900-01-01 00:00:00" -a modulo_end,TIME,c,c,"1901-01-01 00:00:00" $file_out
   ncatted -a calendar,TIME,o,c,"noleap" $file_out
done















