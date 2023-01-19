#!/bin/bash
# this script should run once a day and scheduled with crontab at 0 12 * * *
# clean sd card
path_to_compress=`date '+%Y/%m' --date '2 week ago'``
day_to_compress=`date '+%d' --date '2 week ago'`
path_to_delete=`date '+%Y/%m' --date '4 week ago'`
day_to_delete=`date '+%d' --date '4 week ago'``
#compress images of the previous month
gzip /var/www/html/data/$path_to_compress/*$day_to_compress_*.dng
#remove the old images
rm  -f /var/www/html/data/$path_to_delete/*$day_to_delete_*.dng.gz
rm  -f /var/www/html/data/$path_to_delete/*$day_to_delete_*.jpg
path_to_compress=`date '+%Y/%m' --date '8 week ago'`
day_to_compress=`date '+%d' --date '8 week ago'`
path_to_delete=`date '+%Y/%m' --date '16 week ago'`
day_to_delete=`date '+%d' --date '16 week ago'`
#compress images
gzip /home/sand/data/$path_to_compress/*$day_to_compress_*.dng
#remove the old images
rm  -f /home/sand/data/$path_to_delete/*$day_to_delete_*.dng.gz
rm  -f /home/sand/data/$path_to_delete/*$day_to_delete_*.jpg
# suffix for archiving log file
logsuf=`date '+_%j`
cp -f /home/sand/nightmon.log "nightmon.log_"$logsuf
echo "" > /home/sand/nightmon.log
