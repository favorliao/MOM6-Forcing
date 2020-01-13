#!/bin/bash

ori_dir='/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/data_dfs_mom6'
res_dir='/tigress/GEOCLIM/LRGROUP/datasets/dfs_forcing_drowned/res_dfs_mom6_3yrs'

#declare -a varall=("precip" "q10" "radlw" "radsw" "slp" "snow" "t10" "u10" "v10")
#declare -a varall=("precip" "radlw" "radsw" "slp" "snow" "u10" "v10")
declare -a varall=("snow")
varlength=${#varall[@]}
# here to select the variable for merging
for (( i=1; i<${varlength}+1; i++ ));
#for (( i=1; i<2; i++ ));
do
  var=${varall[$i-1]}
  cd $ori_dir
  fileall=(`ls *$var*.nc`)
  cd $res_dir
  filelength=${#fileall[@]}
  jend=$((${filelength}-1))
  echo "variable " $i "/" ${varlength} " : " $var
# here to select the time period for merging data
#  for (( j=1; j<=${filelength}-1; j++ )); #60-1, to 2016
  for (( j=1; j<=$jend; j++ ));do
#  for (( j=1; j<=2; j++ ));
     if [ "$j" -ge 2 ]; then
        echo $year_end
     #   cdo -O -b f32 mergetime $filename_before $filename_current $filename_next $filename_output
     elif ["$j" -ge $jend]; then
        echo $year_end
     else
        echo $year_end
     #   cdo -O -b f32 mergetime $filename_current $filename_next $filename_output
     fi
    filename_current="$ori_dir/${fileall[$j-1]}"
    filename_next="$ori_dir/${fileall[$j]}"
    filename_output="$res_dir/${fileall[$j-1]}"
    if [ "$j" -ge 2 ]; then 
       filename_before="$ori_dir/${fileall[$j-2]}"
       year_start=$((1957+$j-1))
       echo "filename inbefore "$j "/" ${filelength} " : " $filename_before
    else
       year_start=$((1957+$j))
    fi

    if [ "$j" -eq "$jend" ]; then
       year_end=$((1957+$j+1))
    else
       year_end=$((1957+$j+2))
    fi

    echo "filename incurrent " $j "/" ${filelength} " : " $filename_current
    echo "filename innext "$j "/" ${filelength} " : " $filename_next
    echo "filename out  "$j "/" ${filelength} " : " $filename_output
    echo $year_start
    echo $year_end

     if [ "$j" -ge 2 ]; then 
        echo $year_end
     #   cdo -O -b f32 mergetime $filename_before $filename_current $filename_next $filename_output
     elif ["$j" -ge $jend]; then 
        echo $year_end
     else
        echo $year_end
     #   cdo -O -b f32 mergetime $filename_current $filename_next $filename_output
     fi
  
     echo "     change time att ..."
     #ncatted -O -a modulo,TIME,c,c,' ' -a modulo_beg,TIME,c,c,"$year_start-01-01 00:00:00" -a modulo_end,TIME,c,c,"$year_end-01-01 00:00:00" $filename_output
     echo "     set calendar ..."
     #ncatted  -a  calendar,TIME,o,c,"gregorian" $filename_output
#     cdo -b f32 -f nc copy $filename_output $filename_output
  done

done




