#!/bin/bash
# clean sd card
path_to_compress=`date '+%Y/%m' --date '1 month ago'``
path_to_delete=`date '+%Y/%m' --date '6 month ago'``
#compress images of the previous month
gzip /var/www/html/data/$path_to_compress/*.dng
#remove the old images
rm  -f /var/www/html/data/$path_to_delete/*.dng.gz
