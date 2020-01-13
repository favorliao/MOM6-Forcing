#!/bin/bash

dir='/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/res_dfs_mom6_3yrs'
cd $dir
finall=(`ls dfs_mom6*.nc`)

#dfs_mom6_snow_2005_10oct2018.nc

for fin in "${finall[@]}"
do
   echo "$fin"
   slen=$(echo -n $fin | wc -m)
   var="${fin:0:3}5.2_${fin:4:slen}"
   mv $fin $var
   echo "$var"
done
