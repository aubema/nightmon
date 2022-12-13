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
		 read tv bidon < Current_A_tint.tmp
		 if [ -z "$ta" ]
		 then echo "ta not available, setting it to 1/10s"
		      let ta=100000
		 fi
		 let satmax=1000
		 # capture an image with V camera
		 while [ "$satmax" -gt 99 ] || [ "$satmax" -lt 70 ]
		 do	rm -f capture_1*
		    echo "Ta=" $ta
		 		captureA.py -t $ta -g $gain
				if [ -f "capture_1.dng" ]
				then lisc perc capture_1.dng -p 99.9  > saturation.tmp
				     maxsatpercent.py > capture.tmp
				     read satmax bidon  < capture.tmp
						 echo "satmax=" $satmax
			       if [ "$satmax" -ge 100 ]
			       then  let ta=ta/4
					 elif [ "$satmax" -lt 70 ]
			       then let ta=90*ta/satmax
					   fi
			  else echo "Problem with A camera."
				  	 exit 0
				fi
		 done
		 echo  $ta > Current_A_tint.tmp
		 echo "V integration time: " $ta >> nightmon.log
}

take_pictureB() {
		      #  Take pictures of various integration times starting from a smaller to get the right integration time (max around 0.8)
		 		 echo "Taking B picture"
		 		 read tb bidon < Current_B_tint.tmp
		 		 if [ -z "$tb" ]
		 		 then echo "tb not available, setting it to 1/10s"
		 		      let tb=100000
		 		 fi
		 		 let satmax=1000
		 let satmax=1000
		 while [ "$satmax" -gt 99 ] || [ "$satmax" -lt 70 ]
		 do	rm -f capture_2*
		    echo "Tb=" $tb
		 		captureB.py -t $tb -g $gain
				lisc perc capture_2.dng -p 99.9
				if [ -f "capture_2.dng" ]
		    then lisc perc capture_2.dng -p 99.9  > saturation.tmp
				     maxsatpercent.py > capture.tmp
				     read satmax bidon  < capture.tmp
						 echo "satmax=" $satmax
			       if [ "$satmax" -ge 100 ]
			       then  let tb=tb/4
			       elif [ "$satmax" -lt 70 ]
			       then let tb=90*tb/satmax
					   fi
			  else echo "Problem with V camera."
				  	 exit 0
				fi
		 done
		 echo  $tb > Current_B_tint.tmp
		 echo "R integration time: " $tb >> nightmon.log
}


#
# ==================================
#
# main
#
rm -f *.tmp
user="sand"
gain=8
max_lum=10000  # 1000000 = 1sec
darkimg="dark-gain8-t100000.dng"
# possible bands JV JR R G B (max 2 bands)
bands=(JV JR)
model="RpiHQ-JFilters"  # other choices are "RpiHQ" and "A7S"
cams=(A B)
# Standard extinctions for Observatorio Roque de los Muchachos (0.102 0.0547) (JV JR)
extinct=(0.102 0.0547)
echo ${bands[0]};for i in ${bands[@]};do echo $i;done
echo  "10000 us" > Current_A_tint.tmp
echo  "10000 us" > Current_B_tint.tmp
# get the site name
/bin/grep "SITE" /home/$user/nightmon_config > ligne.tmp
read bidon bidon sitename bidon < ligne.tmp
basepath="/var/www/html/data"
backpath="/home/$user/data"
echo "A shot"
take_pictureA
y=`date +%Y`
mo=`date +%m`
dA=`date +%d`
basenameA=`date +%Y-%m-%d_%H-%M-%S`
baseday=`date +%Y-%m-%d`
basename[0]="$basenameA"
read  tv toto < Current_A_tint.tmp
# writing to logfile
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
echo $y $mo $dA " A " $ta $basepath/$yV/$moV/$basenameA"_A_"$ta"_"$gain".dng" >> $basepath/$y/$mo/nightmon.log
echo "=============================="
# rename pictures
cp -f capture_1.dng $basepath/$y/$mo/$basenameA"_A_"$ta"_"$gain".dng"
cp -f capture_1.jpg $basepath/$y/$mo/$basenameA"_A_"$ta"_"$gain".jpg"
mv capture_1.dng $basenameA"_A_"$ta"_"$gain".dng"

echo "B shot"
take_pictureB
basenameB=`date +%Y-%m-%d_%H-%M-%S`
basename[1]="$basenameB"
read  tb toto < Current_B_tint.tmp
echo $y $mo $d " B " $tb $basenameB"_B_"$tb"_"$gain".dng" >> $basepath/$y/$mo/nightmon.log
echo "=============================="
# rename pictures
cp -f capture_2.dng $basepath/$y/$mo/$basenameB"_B_"$tb"_"$gain".dng"
cp -f capture_2.jpg $basepath/$y/$mo/$basenameB"_B_"$tb"_"$gain".jpg"
mv capture_2.dng $basenameB"_B_"$tb"_"$gain".dng"


# check for the night by reading the latest optimal integration time
if [ $ta -lt $max_lum ]
then echo "Too much light. It is probably daytime."
     move_cams.py 2000 1
		 move_cams.py -1500 1
		 echo "Let's keep the camera inside for 15 min"
		 exit 0
else move_cams.py 2000 1
fi

echo "=============================="
# process sky IMAGES
let n=0
for b in ${bands[@]}
do 	if [ $n -eq 0 ]
    then let t=ta
		else let t=tb
		fi
		python3 /usr/local/bin/ProcessNightMon.py -s ${basename[$n]}"_"${cams[$n]}"_"$t"_"$gain".dng" -d /home/$user/git/nightmon/data/Darks/$darkimg -b $b -e ${extinct[$n]} -c ${cams[$n]} -m $model
		if [ -f $band"calibration"${basename[$n]}".png" ]
		then
			# rename plots
			mv $band"calibration"${basename[$n]}".png" $basepath/$y/$mo/
			mv $band"_calSbBkg_"${basename[$n]}".png" $basepath/$y/$mo/
			mv $band"_calSbTot_"${basename[$n]}".png" $basepath/$y/$mo/
			mv $band"_Stars_Match_"${basename[$n]}".png" $basepath/$y/$mo/

			# backup plots
			cp -f $basepath/$y/$mo/$b"calibration"${basename[$n]}".png" $backpath/$y/$mo/
			cp -f $basepath/$y/$mo/$b"_calSbBkg_"${basename[$n]}".png" $backpath/$y/$mo/
			cp -f $basepath/$y/$mo/$b"_calSbTot_"${basename[$n]}".png" $backpath/$y/$mo/
			cp -f $basepath/$y/$mo/$band"_Stars_Match_"${basename[$n]}".png" $backpath/$y/$mo/
		fi

		# backup output files
		cp -f $basepath/$y/$mo/nightmon.log $backpath/$y/$mo/
		if [ -f $basepath/$y/$mo/"calibrated_"$baseday"_sky.csv" ]
		then cat "calibrated_"$b"_"$baseday"_sky.csv" | grep -v "Loc_Name" | grep -v "(pixel)"  >> $basepath/$y/$mo/"calibrated_"$baseday"_sky.csv"
		else cat "calibrated_"$b"_"$baseday"_sky.csv" | grep "Loc_Name" | grep "(pixel)"  > $basepath/$y/$mo/"calibrated_"$baseday"_sky.csv"
		fi
		cp -f $basepath/$y/$mo/"calibrated_"$baseday"_sky.csv" $backpath/$y/$mo/

		# clean directory
		rm ${basename[$n]}"_"${cams[$n]}"_"$t"_"$gain".dng"
		rm "calibrated_"$b"_"$baseday"_sky.csv"
		let n=n+1
done
