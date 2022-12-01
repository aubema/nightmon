#!/bin/bash
# clean sd card
path_to_clean=`date '+%Y/%m' --date '1 month ago'``
#remove the oldest images
rm  ls -l --sort=time /var/www/html/data/$path_to_clean/*_R_*.dng | sed -n 2p  | awk '{print $NF}'
rm  ls -l --sort=time /var/www/html/data/$path_to_clean/*_V_*.dng | sed -n 2p  | awk '{print $NF}'
rm  ls -l --sort=time /var/www/html/data/$path_to_clean/*_V_*.jpg | sed -n 2p  | awk '{print $NF}'
rm  ls -l --sort=time /var/www/html/data/$path_to_clean/*_V_*.jpg | sed -n 2p  | awk '{print $NF}'
