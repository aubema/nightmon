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
take_pictureV() {
     #  Take pictures of various integration times starting from a smaller to get the right integration time (max around 0.8)
		 echo "Taking V picture"
		 read tv bidon < Current_V_tint.tmp
		 if [ -z "$tv" ]
		 then echo "tv not available, setting it to 1/10s"
		      let tv=100000
		 fi
		 let satmax=1000
		 # capture an image with V camera
		 while [ "$satmax" -gt 90 ] || [ "$satmax" -lt 70 ]
		 do	rm -f capture_1*
		    echo "Tv=" $tv
		 		captureA.py -t $tv -g $gain
				if [ -f "capture_1.dng" ]
				then lisc perc capture_1.dng -p 99.9  > saturation.tmp
				     maxsatpercent.py > capture.tmp
				     read satmax bidon  < capture.tmp
						 echo "satmax=" $satmax
			       if [ "$satmax" -ge 100 ]
			       then  let tv=tv/10
					 elif [ "$satmax" -lt 70 ]
			       then let tv=90*tv/satmax
					   fi
			  else echo "Problem with V camera."
				  	 exit 0
				fi
		 done
		 echo  $tv > Current_V_tint.tmp
		 echo "V integration time: " $tv >> nightmon.log
}

take_pictureR() {
		      #  Take pictures of various integration times starting from a smaller to get the right integration time (max around 0.8)
		 		 echo "Taking pictures"
		 		 read tr bidon < Current_R_tint.tmp
		 		 if [ -z "$tr" ]
		 		 then echo "tr not available, setting it to 1/10s"
		 		      let tr=100000
		 		 fi
		 		 let satmax=1000
		 let satmax=1000
		 while [ "$satmax" -gt 90 ] || [ "$satmax" -lt 70 ]
		 do	rm -f capture_2*
		    echo "Tr=" $tr
		 		captureB.py -t $tr -g $gain
				lisc perc capture_2.dng -p 99.9
				if [ -f "capture_2.dng" ]
		    then lisc perc capture_2.dng -p 99.9  > saturation.tmp
				     maxsatpercent.py > capture.tmp
				     read satmax bidon  < capture.tmp
						 echo "satmax=" $satmax
			       if [ "$satmax" -ge 100 ]
			       then  let tr=tr/10
			       elif [ "$satmax" -lt 70 ]
			       then let tr=90*tr/satmax
					   fi
			  else echo "Problem with V camera."
				  	 exit 0
				fi
		 done
		 echo  $tr > Current_R_tint.tmp
		 echo "R integration time: " $tr >> nightmon.log
}


#
# ==================================
#
# main
#
gain=8
max_lum=10000  # 1000000 = 1sec
darkimg="dark-gain8-t100000.dng"
echo  "10000 us" > Current_V_tint.tmp
echo  "10000 us" > Current_R_tint.tmp
# get the site name
/bin/grep "SITE" /home/sand/nightmon_config > ligne.tmp
read bidon bidon sitename bidon < ligne.tmp

basepath="/var/www/html/data"
backpath="/home/sand/data"
echo "Johnson V shot"
take_pictureV
yV=`date +%Y`
moV=`date +%m`
dV=`date +%d`
basenameV=`date +%Y-%m-%d_%H-%M-%S`
basedayV=`date +%Y-%m-%d`
read  tv toto < Current_V_tint.tmp
# writing to logfile
if [ ! -d $basepath/$yV ]
then mkdir $basepath/$yV
fi
if [ ! -d $basepath/$yV/$moV ]
then /bin/mkdir $basepath/$yV/$moV
fi
echo $yV $moV $dV $h $mi $s " V " $tv $basepath/$yV/$moV/$basenameV"_V_"$tv"_"$gain".dng" >> $basepath/$yV/$moV/nightmon.log
echo "=============================="
# rename pictures
cp -f capture_1.dng $basepath/$yV/$moV/$basenameV"_V_"$tv"_"$gain".dng"
cp -f capture_1.jpg $basepath/$yV/$moV/$basenameV"_V_"$tv"_"$gain".jpg"
mv capture_1.dng $basenameV"_V_"$tv"_"$gain".dng"

echo "Johnson R shot"
take_pictureR
yR=`date +%Y`
moR=`date +%m`
dR=`date +%d`
basenameR=`date +%Y-%m-%d_%H-%M-%S`
read  tr toto < Current_R_tint.tmp
# writing to logfile
if [ ! -d $basepath/$yR ]
then mkdir $basepath/$yR
fi
if [ ! -d $basepath/$yR/$moR ]
then /bin/mkdir $basepath/$yR/$moR
fi
echo $yR $moR $dR $h $mi $s " R " $tr $basenameR"_R_"$tv"_"$gain".dng" >> $basepath/$yR/$moR/nightmon.log
echo "=============================="
# rename pictures
cp -f capture_2.dng $basepath/$yR/$moR/$basenameR"_R_"$tr"_"$gain".dng"
cp -f capture_2.jpg $basepath/$yR/$moR/$basenameR"_R_"$tr"_"$gain".jpg"
mv capture_2.dng $basenameR"_R_"$tr"_"$gain".dng"


# check for the night by reading the latest optimal integration time
if [ $tv -lt $max_lum ]
then echo "Too much light. It is probably daytime."
     move_cams.py 2000 1
		 move_cams.py -1500 1
		 echo "Let's keep the camera inside for 15 min"
else move_cams.py 2000 1
fi

echo "=============================="
# process sky IMAGES
echo
python3 /usr/local/bin/ProcessNightMon-JVR2H.py -v $basenameV"_V_"$tv"_"$gain".dng" -r $basenameR"_R_"$tr"_"$gain".dng" -d /home/sand/git/nightmon/data/Darks/$darkimg
# rename pictures
mv Vzeropoint_corr$basenameV.png $basepath/$yV/$moV/
mv Rzeropoint_corr$basenameV.png $basepath/$yR/$moR/Rzeropoint_corr$basenameR.png
mv VcalSbBkg$basenameV.png $basepath/$yV/$moV/
mv RcalSbBkg$basenameV.png $basepath/$yR/$moR/RcalSbBkg$basenameR.png
mv VStars_Match$basenameV.png $basepath/$yV/$moV/
mv RStars_Match$basenameV.png $basepath/$yR/$moR/RStars_Match$basenameR.png
cat "calibrated_"$basedayV"_sky.csv" | grep -v "Loc_Name" | grep -v "(pixel)"  >> $basepath/$yV/$moV/"calibrated_"$basedayV"_sky.csv"
# backup important files
if [ ! -d $backpath/$yV ]
then mkdir $backpath/$yV
fi
if [ ! -d $backpath/$yV/$moV ]
then /bin/mkdir $backpath/$yV/$moV
fi
if [ ! -d $backpath/$yR ]
then mkdir $backpath/$yR
fi
if [ ! -d $backpath/$yR/$moR ]
then /bin/mkdir $backpath/$yR/$moR
fi
rm $basenameR"_R_"$tr"_"$gain".dng"
rm $basenameR"_V_"$tv"_"$gain".dng"
cp -f $basepath/$yV/$moV/Vzeropoint_corr$basenameV.png $backpath/$yV/$moV/
cp -f $basepath/$yR/$moR/Rzeropoint_corr$basenameR.png $backpath/$yR/$moR/
cp -f $basepath/$yV/$moV/VStars_Match$basenameV.png $backpath/$yV/$moV/
cp -f $basepath/$yR/$moR/RStars_Match$basenameR.png $backpath/$yR/$moR/
cp -f $basepath/$yV/$moV/nightmon.log $backpath/$yV/$moV/
cp -f $basepath/$yR/$moR/nightmon.log $backpath/$yR/$moR/
cp -f $basepath/$yV/$moV/"calibrated_"$basedayV"_sky.csv" $backpath/$yV/$moV/
