#!/bin/bash
#
#
#    Copyright (C) 2022  Martin Aube
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Contact: martin.aube@cegepsherbrooke.qc.ca
#

#
# =============================
# find appropriate integration time
#
take_pictures() {
	nightt=120000000   # 2 minutes
	nightg=16
	#  Take pictures of various integration times starting from a smaller to get the right integration time (max around 0.8)
	echo "Taking A picture"
	echo $path
	let ta=nightt
	let gain=nightg
	let satmax=1000
	rm -f $path"/capture_1*"
	rm -f $path"/capture_3*"
	echo "Shooting "$ta" micro seconds... with gain " $gain 
	/usr/bin/python3 /usr/local/bin/captureA.py -t $ta -g $gain
	if [ -f $path"/capture_1.dng" ] ; then
		/usr/local/bin/lisc perc $path"/capture_1.dng" -p 99.9  > $path"/saturation.tmp"
		/usr/bin/python3 /usr/local/bin/maxsatpercent.py > $path"/capture.tmp"
		read satmax bidon  < $path"/capture.tmp"
		echo "satmax=" $satmax
		while [ "$satmax" -ge 80 ] && [ "$ta" -gt 1200 ]
		do let ta=ta/10
			rm -f $path"/capture_1*"
			echo "Shooting "$ta" micro seconds... with gain " $gain
			/usr/bin/python3 /usr/local/bin/captureA.py -t $ta -g $gain
			if [ -f $path"/capture_1.dng" ] ; then
				/usr/local/bin/lisc perc $path"/capture_1.dng" -p 99.9  > $path"/saturation.tmp"
				/usr/bin/python3 /usr/local/bin/maxsatpercent.py > $path"/capture.tmp"
				read satmax bidon  < $path"/capture.tmp"
				echo "satmax=" $satmax
			else
				echo "Problem with camera A."
				exit 0
			fi
		done
		if [ "$satmax" -ge 80 ]
		then let gain=2
			  rm -f $path"/capture_1*"
			  echo "Shooting "$ta" micro seconds... with gain " $gain
			  /usr/bin/python3 /usr/local/bin/captureA.py -t $ta -g $gain
			  if [ -f $path"/capture_1.dng" ] ; then
				  /usr/local/bin/lisc perc $path"/capture_1.dng" -p 99.9  > $path"/saturation.tmp"
				  /usr/bin/python3 /usr/local/bin/maxsatpercent.py > $path"/capture.tmp"
				  read satmax bidon  < $path"/capture.tmp"
				  echo "satmax=" $satmax
			  else
			     echo "Problem with camera A."
				  exit 0
			  fi
		fi
		/usr/bin/python3 /usr/local/bin/captureC.py -t $ta -g $gain
	else
		echo "Problem with camera A."
		exit 0
	fi
	echo  $ta > $path"/Current_tint.tmp"
	echo  $gain > $path"/Current_gain.tmp"
	# flush ram cache to correct a memory leak in the camera library
	/usr/bin/sync
	/usr/bin/echo 3 > /proc/sys/vm/drop_caches	
}

# ==================================
# global positioning system
globalpos () {

     rm -f /root/*.tmp
     bash -c '/usr/bin/gpspipe -w -n 2 | sed -e "s/,/\n/g" | grep activated | tail -1 | sed "s/n\"/ /g" |sed -e "s/\"/ /g" |  sed -e"s/activated//g" | sed -e "s/ //g" > /home/sand/coords.tmp'
     read gpstime1 < /home/sand/coords.tmp
     gpstime1="${gpstime1:1}"
     bash -c '/usr/bin/gpspipe -w -n 4 | sed -e "s/,/\n/g" | grep time | tail -1 | sed "s/n\"/ /g" |sed -e "s/\"/ /g" |  sed -e"s/time//g" | sed -e "s/ //g" > /home/sand/coords.tmp'
     read gpstime2 < /home/sand/coords.tmp
     gpstime2="${gpstime2:1}"
     echo "t" $gpstime1 $gpstime2
     sec1=`/usr/bin/date -d "$gpstime1" +%s`
     sec2=`/usr/bin/date -d "$gpstime2" +%s`
     if [ $sec1 -gt $sec2 ]
     then gpstime=$gpstime1
     else gpstime=$gpstime2
     fi
     echo $gpstime
     bash -c '/usr/bin/gpspipe -w -n 5 | sed -e "s/,/\n/g" | grep lat | tail -1 | sed "s/n\"/ /g" |sed -e "s/\"/ /g" | sed -e "s/:/ /g" | sed -e"s/lat//g" | sed -e "s/ //g" > /home/sand/coords.tmp'
     read lat < /home/sand/coords.tmp
     bash -c '/usr/bin/gpspipe -w -n 5 | sed -e "s/,/\n/g" | grep lon | tail -1 | sed "s/n\"/ /g" |sed -e "s/\"/ /g" | sed -e "s/:/ /g" | sed -e "s/lo//g" | sed -e "s/ //g" > /home/sand/coords.tmp'
     read lon < /home/sand/coords.tmp
     bash -c '/usr/bin/gpspipe -w -n 5 | sed -e "s/,/\n/g" | grep alt | tail -1 | sed "s/n\"/ /g" |sed -e "s/\"/ /g" | sed -e "s/:/ /g" | sed -e "s/alt//g" | sed -e "s/ //g" > /home/sand/coords.tmp'
     read alt < /home/sand/coords.tmp
     echo $lat $lon $alt
     if [ -z "${lon}" ]
     then let lon=0
          let lat=0
          let alt=0
     fi 
     # /bin/echo "GPS gives Latitude:" $lat ", Longitude:" $lon "and Altitude:" $alt
     /bin/echo "Lat.:" $lat ", Lon.:" $lon " Alt.:" $alt  > /home/sand/gps.log
     echo $gpstime > /home/sand/date_gps.log
}


#
# ==================================
#
# main
#
# setting constants
user="sand"
zpoint=1.0
gain=16
min_gain=5  # if below the image will not be calibrated
darkimg="dark-gain16-t100000.dng"
bands=(JB JV JR)
model="RpiHQ"  # other choices are "RpiHQ" and "A7S"
cams=(A C)
basepath="/var/www/html/data"
backpath="/home/"$user"/data"
path="/home/"$user
# get the site name
/bin/grep "SITE" $path"/nightmon_config" > $path"/ligne.tmp"
read bidon bidon sitename bidon < $path"/ligne.tmp"
# wait 2 min to start (enough time for ntp sync)
# echo "Waiting 2 min before starting measurements..."
# /bin/sleep 120
#OPTIONS
gopt=0
while getopts 'k:' OPTION ; do
	case $OPTION in
	k) kopt=1
	Kalmode=$OPTARG
	;;
	esac
done
# Reading calibration mode (stars or fixed)
if [ "$gopt" != "1" ] ; then
	Kalmode=$1
fi



while : ; do
	processflag=0
	time1=`date +%s`
	# globalpos
	# secg=`/usr/bin/date -d "$gpstime" +%s`
	# if [ $secg -gt $time1 ] ; then
	#	/usr/bin/date -s $gpstime
	# fi
	echo "Shooting..."
	take_pictures
	y=`date --date="2 minutes ago" +%Y`
	mo=`date --date="2 minutes ago" +%m`
	dA=`date --date="2 minutes ago" +%d`
	basenameA=`date --date="2 minutes ago" +%Y-%m-%d_%H-%M-%S`
	basenameC=`date +%Y-%m-%d_%H-%M-%S`
	echo $basenameA
	baseday=`date --date="2 minutes ago" +%Y-%m-%d`
	basename[0]="$basenameA"
	# find the closest extinction in time from /home/sand/extinction_data
	Emindelay=10000000000
	while read -r line; do
		echo $line > $path"/ligne.tmp"
		read Ye Me De Ve Re bidon < $path"/ligne.tmp"
		if [ $Ye != "#" ] ; then
			timee=`date --date=$Ye"-"$Me"-"$De" 00:00:01" +%s`
			let DTe=time1-timee
			DTe="${DTe/#-}"   # absolute value
			if [ $DTe -lt $Emindelay ] ; then
				let Emindelay=DTe
				extinct[0]=$Ve
				extinct[1]=$Re
			fi
		fi
	done < $path"/extinction_data"

	read  tv toto < $path"/Current_tint.tmp"
	read  gain toto < $path"/Current_gain.tmp"
	if [ ! -d $basepath/$y ] ; then
		mkdir $basepath/$y
	fi
	if [ ! -d $basepath/$y/$mo ] ; then
		/bin/mkdir $basepath/$y/$mo
	fi
	if [ ! -d $backpath/$y ] ; then
		mkdir $backpath/$y
	fi
	if [ ! -d $backpath/$y/$mo ] ; then
		/bin/mkdir $backpath/$y/$mo
	fi
	echo "=============================="
	# rename pictures
	if [ $gain -ge $min_gain ] ; then
		# daytime only keep the jpg
		cp -f $path"/capture_1.dng" $basepath/$y/$mo/$basenameA"_A_"$ta"_"$gain".dng"
		cp -f $path"/capture_1.dng" $backpath/$y/$mo/$basenameA"_A_"$ta"_"$gain".dng"
		cp -f $path"/capture_1.dng" $path/$basenameA"_A_"$ta"_"$gain".dng"
		processflag=1
	fi
	cp -f $path"/capture_1.jpg" $basepath/$y/$mo/$basenameA"_A_"$ta"_"$gain".jpg"
	cp -f $path"/capture_1.jpg" $backpath/$y/$mo/$basenameA"_A_"$ta"_"$gain".jpg"
	cp -f $path"/capture_1.jpg" $path/$basenameA"_A_"$ta"_"$gain".jpg"
	basename[1]="$basenameC"
	read  tb toto < $path"/Current_tint.tmp"
	echo "=============================="
	# rename pictures
	if [ $ta -ge $max_lum ] ; then
		cp -f $path"/capture_3.dng" $basepath/$y/$mo/$basenameC"_C_"$tb"_"$gain".dng"
		cp -f $path"/capture_3.dng" $backpath/$y/$mo/$basenameC"_C_"$tb"_"$gain".dng"
		cp -f $path"/capture_3.dng" $path/$basenameC"_C_"$tb"_"$gain".dng"
	fi
	cp -f $path"/capture_3.jpg" $basepath/$y/$mo/$basenameC"_C_"$tb"_"$gain".jpg"
	cp -f $path"/capture_3.jpg" $backpath/$y/$mo/$basenameC"_C_"$tb"_"$gain".jpg"
	cp -f $path"/capture_3.jpg" $path/$basenameC"_C_"$tb"_"$gain".jpg"
	echo "=============================="
	if [ $processflag -eq 1 ] ; then
		# process sky IMAGES
		let n=0
		for p in ${cams[@]} ; do 
			for b in ${bands[@]} ; do	
			
				if [ $n -eq 0 ] ; then
				let t=ta
				else
					let t=tb
        		fi
				# determine the zeropoint according to the integration time and camera
				#
				#
				#
				/usr/bin/python3 /usr/local/bin/ProcessNightMon.py -i $basepath/$y/$mo/${basename[$n]}"_"${cams[$n]}"_"$t"_"$gain".dng" -d $path"/git/nightmon/data/Darks/"$darkimg -b $b -k fixed -m $model
				if [ -f $path"/"$p"_"$b"_calibration_"${basename[$n]}".png" ] ; then
					mv $path"/"$p"_"$b"_calibration_"${basename[$n]}".png" $basepath/$y/$mo/
					cp -f $basepath"/"$y"/"$mo"/"$p"_"$b"_calibration_"${basename[$n]}".png" $backpath"/"$y"/"$mo"/"
				fi
				if [ -f $path"/"$p"_"$b"_calSbBkg_"${basename[$n]}".png" ] ; then
					mv $path"/"$p"_"$b"_calSbBkg_"${basename[$n]}".png" $basepath/$y/$mo/
					cp -f $basepath"/"$y"/"$mo"/"$p"_"$b"_calSbBkg_"${basename[$n]}".png" $backpath"/"$y"/"$mo"/"
				fi
				if [ -f $path"/"$p"_"$b"_calSbTot_"${basename[$n]}".png" ] ; then
					mv $path"/"$p"_"$b"_calSbTot_"${basename[$n]}".png" $basepath/$y/$mo/
					cp -f $basepath"/"$y"/"$mo"/"$p"_"$b"_calSbTot_"${basename[$n]}".png" $backpath"/"$y"/"$mo"/"
				fi
				if [ -f $path"/"$p"_"$b"_Stars_Match_"${basename[$n]}".png" ] ; then
					mv $path"/"$p"_"$b"_Stars_Match_"${basename[$n]}".png" $basepath/$y/$mo/
					cp -f $basepath"/"$y"/"$mo"/"$p"_"$b"_Stars_Match_"${basename[$n]}".png" $backpath"/"$y"/"$mo"/"
				fi
				# backup calibration file
				if [ -f $path"/"$p"_"$b"_calibration_stars_"${basename[$n]}".csv" ] ; then 
					mv $path"/"$p"_"$b"_calibration_stars_"${basename[$n]}".csv" $basepath/$y/$mo/
					cp -f $basepath"/"$y"/"$mo"/"$p"_"$b"_calibration_stars_"${basename[$n]}".csv" $backpath"/"$y"/"$mo"/"
				fi




				# backup output files
				if [ -f $basepath"/"$y"/"$mo"/calibrated_"$baseday"_sky.csv" ] ; then
					cat $path"/calibrated_"$p"_"$b"_"$baseday"_sky.csv" | grep -v "Loc_Name" | grep -v "(pixel)"  >> $basepath"/"$y"/"$mo"/calibrated_"$baseday"_sky.csv"
         	else
        		 			cat $path"/calibrated_"$p"_"$b"_"$baseday"_sky.csv"  > $basepath"/"$y"/"$mo"/calibrated_"$baseday"_sky.csv"
				fi
				cp -f $basepath"/"$y"/"$mo"/calibrated_"$baseday"_sky.csv" $backpath"/"$y"/"$mo"/"
			done
			let n=n+1
			# clean directory
			rm $path"/"${basename[$n]}"_"${cams[$n]}"_"$t"_"$gain".dng"
			rm $path"/calibrated_"$p"_"$b"_"$baseday"_sky.csv"
		done
	fi
	time2=`date +%s`
	let idle=900-time2+time1  # one measurement every 15 min (15*60=900)
	if [ $idle -lt 0 ] ; then 
		let idle=0
	fi
	echo "Wait " $idle "s before next reading."
	/bin/sleep $idle
	time1=`date +%s`
	rm -f $path/*.dng
	rm -f $path/*.jpg
done
