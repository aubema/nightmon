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
take_pictureA() {
     #  Take pictures of various integration times starting from a smaller to get the right integration time (max around 0.8)
		 echo "Taking A picture"
		 echo $path
		 read ta bidon < $path"/Current_A_tint.tmp"
		 if [ -z "$ta" ] || [ "$ta" -lt 800 ]
		 then echo "ta not available, setting it to 1/10s"
		      let ta=3750000
		 fi
		 let satmax=1000
		 # capture an image with V camera
		 while [ "$satmax" -gt 99 ] || [ "$satmax" -lt 50 ]
		 do	rm -f $path"/capture_1*"
		    echo "Ta=" $ta
		 		/usr/bin/python3 /usr/local/bin/captureA.py -t $ta -g $gain
				if [ "$ta" -eq 120000000 ]
				then break
				fi
				if [ -f $path"/capture_1.dng" ]
				then /usr/local/bin/lisc perc $path"/capture_1.dng" -p 99.9  > $path"/saturation.tmp"
				     /usr/bin/python3 /usr/local/bin/maxsatpercent.py > $path"/capture.tmp"
				     read satmax bidon  < $path"/capture.tmp"
						 echo "satmax=" $satmax
			       if [ "$satmax" -ge 100 ]
			       then  let ta=ta/4
					 elif [ "$satmax" -lt 50 ]
			       then let ta=ta*2
					   fi
						 #if [ "$ta" -gt 120000000 ]
						 #then let ta=120000000
						 #     echo "Ta=" $ta
             #     /usr/bin/python3 /usr/local/bin/captureA.py -t $ta -g $gain
						 #fi
						 if [ "$ta" -le 1000 ]
						 then break
						 fi
			  else echo "Problem with A camera."
				  	 exit 0
				fi
		 done
		 echo  $ta > $path"/Current_A_tint.tmp"
}

take_pictureB() {
		      #  Take pictures of various integration times starting from a smaller to get the right integration time (max around 0.8)
		 		 echo "Taking B picture"
		 		 read tb bidon < $path"/Current_B_tint.tmp"
		 		 if [ -z "$tb" ]  || [ "$tb" -lt 800 ]
		 		 then echo "tb not available, setting it to 1/10s"
		 		      let tb=3750000
		 		 fi
		 		 let satmax=1000
		 let satmax=1000
		 while [ "$satmax" -gt 99 ] || [ "$satmax" -lt 50 ]
		 do	rm -f $path"/capture_2*"
		    echo "Tb=" $tb
		 		/usr/bin/python3 /usr/local/bin/captureB.py -t $tb -g $gain
				if [ "$tb" -eq 120000000 ]
	 		  then break
			  fi
				/usr/local/bin/lisc perc $path"/capture_2.dng" -p 99.9
				if [ -f $path"/capture_2.dng" ]
		    then /usr/local/bin/lisc perc $path"/capture_2.dng" -p 99.9  > $path"/saturation.tmp"
				     /usr/bin/python3 /usr/local/bin/maxsatpercent.py > $path"/capture.tmp"
				     read satmax bidon  < $path"/capture.tmp"
						 echo "satmax=" $satmax
			       if [ "$satmax" -ge 100 ]
			       then  let tb=tb/4
					 elif [ "$satmax" -lt 50 ]
			       then let tb=tb*2
					   fi
						 #if [ "$tb" -gt 120000000 ]
						 #then let tb=120000000
						 #     echo "Tb=" $tb
             #     /usr/bin/python3 /usr/local/bin/captureB.py -t $tb -g $gain
						 #fi
						 if [ "$tb" -le 1000 ]
						 then break
					   fi
			  else echo "Problem with B camera."
				  	 exit 0
				fi
		 done
		 echo  $tb > 		$path"/Current_B_tint.tmp"
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
max_lum=100000  # 1000000 = 1sec
darkimg="dark-gain16-t100000.dng"
# possible bands JV JR R G B (max 2 bands)
bands=(JV JR)
model="RpiHQ-JFilters"  # other choices are "RpiHQ" and "A7S"
cams=(A B)
# listen to gpio 5 (limit switch)  state=1 means cams inside, state=0 means cams out
echo "5" > /sys/class/gpio/export
# Standard extinctions for Observatorio Roque de los Muchachos (0.102 0.0547) (JV JR)
extinct=(0.102 0.0547)
basepath="/var/www/html/data"
backpath="/home/"$user"/data"
path="/home/"$user
echo  "3750000 us" > $path"/Current_A_tint.tmp"
echo  "3750000 us" > $path"/Current_B_tint.tmp"
# get the site name
/bin/grep "SITE" $path"/nightmon_config" > $path"/ligne.tmp"
read bidon bidon sitename bidon < $path"/ligne.tmp"
# wait 2 min to start (enough time for ntp sync)
echo "Waiting 2 min before starting measurements..."
/bin/sleep 120
#OPTIONS
gopt=0
while getopts 'k:' OPTION
do
   case $OPTION in
      k) kopt=1
         Kalmode=$OPTARG
      ;;
   esac
done
# Reading calibration mode (stars or fixed)
if [ "$gopt" != "1" ]
   then Kalmode=$1
fi



while :
do time1=`date +%s`
   echo "A shot"
   take_pictureA
   y=`date +%Y`
   mo=`date +%m`
   dA=`date +%d`
   basenameA=`date +%Y-%m-%d_%H-%M-%S`
   echo $basenameA
   baseday=`date +%Y-%m-%d`
   basename[0]="$basenameA"
   read  tv toto < $path"/Current_A_tint.tmp"
   if [ ! -d $basepath/$y ]
   then mkdir $basepath/$y
   fi
   if [ ! -d $basepath/$y/$mo ]
   then /bin/mkdir $basepath/$y/$mo
   fi
   if [ ! -d $backpath/$y ]
   then mkdir $backpath/$y
   fi
   if [ ! -d $backpath/$y/$mo ]
   then /bin/mkdir $backpath/$y/$mo
   fi
   echo "=============================="
   # rename pictures
	 if [ $ta -ge $max_lum ]
	 then
      cp -f $path"/capture_1.dng" $basepath/$y/$mo/$basenameA"_A_"$ta"_"$gain".dng"
	    cp -f $path"/capture_1.dng" $backpath/$y/$mo/$basenameA"_A_"$ta"_"$gain".dng"
	    cp -f $path"/capture_1.dng" $path/$basenameA"_A_"$ta"_"$gain".dng"
   fi
   cp -f $path"/capture_1.jpg" $basepath/$y/$mo/$basenameA"_A_"$ta"_"$gain".jpg"
   cp -f $path"/capture_1.jpg" $backpath/$y/$mo/$basenameA"_A_"$ta"_"$gain".jpg"
	 cp -f $path"/capture_1.jpg" $path/$basenameA"_A_"$ta"_"$gain".jpg"
   echo "B shot"
   take_pictureB
   basenameB=`date +%Y-%m-%d_%H-%M-%S`
   basename[1]="$basenameB"
   read  tb toto < $path"/Current_B_tint.tmp"
   echo "=============================="
   # rename pictures
	 if [ $ta -ge $max_lum ]
	 then
      cp -f $path"/capture_2.dng" $basepath/$y/$mo/$basenameB"_B_"$tb"_"$gain".dng"
	    cp -f $path"/capture_2.dng" $backpath/$y/$mo/$basenameB"_B_"$tb"_"$gain".dng"
	    cp -f $path"/capture_2.dng" $path/$basenameB"_B_"$tb"_"$gain".dng"
   fi
   cp -f $path"/capture_2.jpg" $basepath/$y/$mo/$basenameB"_B_"$tb"_"$gain".jpg"
   cp -f $path"/capture_2.jpg" $backpath/$y/$mo/$basenameB"_B_"$tb"_"$gain".jpg"
	 cp -f $path"/capture_2.jpg" $path/$basenameB"_B_"$tb"_"$gain".jpg"
   # check for the night by reading the latest optimal integration time
	 limit=`cat /sys/class/gpio/gpio5/value`
   if [ $ta -lt $max_lum ]
   then echo "Too much light. It is probably daytime."
	    if [ "$limit" == "0" ]
			then
				 echo "Moving cameras..."
         /usr/bin/python3 /usr/local/bin/move_cams.py -1400 1
			fi
      echo "Let's keep the camera inside for 15 min"
			processflag=0
   else /usr/bin/python3 /usr/local/bin/move_cams.py 2000 1
		  processflag=1
   fi
   echo "=============================="
	 if [ $processflag -eq 1 ]
	 then
      # process sky IMAGES
      let n=0
      for b in ${bands[@]}
      do if [ $n -eq 0 ]
         then let t=ta
         else let t=tb
         fi
				 # determine the zeropoint according to the integration time and camera
				 #
				 #
				 #

         /usr/bin/python3 /usr/local/bin/ProcessNightMon.py -i ${basename[$n]}"_"${cams[$n]}"_"$t"_"$gain".dng" -d $path"/git/nightmon/data/Darks/"$darkimg -b $b -e ${extinct[$n]} -c ${cams[$n]} -m $model  -k stars -z $zpoint
         if [ -f $path"/"$b"_calibration_"${basename[$n]}".png" ]
         then
				    mv $path"/"$b"_calibration_"${basename[$n]}".png" $basepath/$y/$mo/
					  cp -f $basepath"/"$y"/"$mo"/"$b"_calibration_"${basename[$n]}".png" $backpath"/"$y"/"$mo"/"
				 fi
				 if [ -f $path"/"$b"_calSbBkg_"${basename[$n]}".png" ]
				 then
					 mv $path"/"$b"_calSbBkg_"${basename[$n]}".png" $basepath/$y/$mo/
					 cp -f $basepath"/"$y"/"$mo"/"$b"_calSbBkg_"${basename[$n]}".png" $backpath"/"$y"/"$mo"/"
			   fi
				 if [ -f $path"/"$b"_calSbTot_"${basename[$n]}".png" ]
				 then
					 mv $path"/"$b"_calSbTot_"${basename[$n]}".png" $basepath/$y/$mo/
					 cp -f $basepath"/"$y"/"$mo"/"$b"_calSbTot_"${basename[$n]}".png" $backpath"/"$y"/"$mo"/"
			   fi
				 if [ -f $path"/"$b"_Stars_Match_"${basename[$n]}".png" ]
				 then
					 mv $path"/"$b"_Stars_Match_"${basename[$n]}".png" $basepath/$y/$mo/
					 cp -f $basepath"/"$y"/"$mo"/"$b"_Stars_Match_"${basename[$n]}".png" $backpath"/"$y"/"$mo"/"
			   fi

         # backup output files
         if [ -f $basepath"/"$y"/"$mo"/calibrated_"$baseday"_sky.csv" ]
         then cat $path"/calibrated_"$b"_"$baseday"_sky.csv" | grep -v "Loc_Name" | grep -v "(pixel)"  >> $basepath"/"$y"/"$mo"/calibrated_"$baseday"_sky.csv"
         else cat $path"/calibrated_"$b"_"$baseday"_sky.csv"  > $basepath"/"$y"/"$mo"/calibrated_"$baseday"_sky.csv"
         fi
         cp -f $basepath"/"$y"/"$mo"/calibrated_"$baseday"_sky.csv" $backpath"/"$y"/"$mo"/"
         # clean directory
         rm $path"/"${basename[$n]}"_"${cams[$n]}"_"$t"_"$gain".dng"
         rm $path"/calibrated_"$b"_"$baseday"_sky.csv"
         let n=n+1
      done
   fi
   time2=`date +%s`
   let idle=900-time2+time1  # one measurement every 15 min (15*60=900)
   if [ $idle -lt 0 ] ; then let idle=0; fi
   echo "Wait " $idle "s before next reading."
   /bin/sleep $idle
	 time1=`date +%s`
	 rm -f $path/*.dng
	 rm -f $path/*.jpg
done
