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
		 if [ -z "$tr" ]
		 then echo "tr not available, setting it to 1/10s"
		      let tr=100000
		 fi
		 let satmax=100
		 # capture an image with V camera
		 while [ $satmax -gt 80 ]
		 do	rm -f capture_1*
		 		captureA.py -t $tv -g $gain
				if [ -f "capture_1.dng" ]
		    then lisc perc capture_1.dng -p 99.9 | sed -e 's/=//g' | sed -e 's/R//g' | sed -e 's/G//g' | sed -e 's/B//g'| sed -e 's/\./ /g' > capture.tmp
		 	   	 read unitr decr unitg decg unitb decb bidon  < capture.tmp
     		   # remove leading zero to the sky brightness
     		   if [ ${decr:0:1} == 0 ]
     		   then decr=`echo $decr | sed -e 's/0//g'`
     		   fi
     		   if [ ${decg:0:1} == 0 ]
     		   then decg=`echo $decg | sed -e 's/0//g'`
     		   fi
     	     if [ ${decb:0:1} == 0 ]
     		   then decb=`echo $decb | sed -e 's/0//g'`
     		   fi
     		   let satr=unitr*1000+decr
	   		   let satg=unitg*1000+decg
     		   let satb=unitb*1000+decb
		 		   sat=($satr $satg $satg)
		 	   	 IFS=$'\n'
		 		   satmax = `echo "${sat[*]}" | sort -nr | head -n1`
		 		   if (satmax -ge 1000)
		       then let tv=tv/2
				   elif (satmax -lt 700)
			     then let tv=80*tv/satmax
			     fi
			  else echo "Problem with V camera."
					 exit 0
				fi
		 done
		 echo  $tv > Current_V_tint.tmp
		 echo "V integration time: " $tv >> nightmon.log
}

take_pictureVR() {
		      #  Take pictures of various integration times starting from a smaller to get the right integration time (max around 0.8)
		 		 echo "Taking pictures"
		 		 read tr bidon < Current_R_tint.tmp
		 		 if [ -z "$tr" ]
		 		 then echo "tr not available, setting it to 1/10s"
		 		      let tr=100000
		 		 fi
		 		 let satmax=100
		 let satmax=100
		 while [ $satmax -gt 80 ]
		 do	rm -f capture_2*
		 		captureB.py -t $tr -g $gain
				if [ -f "capture_2.dng" ]
				then lisc perc capture_2.dng -p 99.9 | sed -e 's/=//g' | sed -e 's/R//g' | sed -e 's/G//g' | sed -e 's/B//g'| sed -e 's/\./ /g' > capture.tmp
				 		read unitr decr unitg decg unitb decb bidon  < capture.tmp
						# remove leading zero to the sky brightness
      		   if [ ${decr:0:1} == 0 ]
      		   then decr=`echo $decr | sed -e 's/0//g'`
      		   fi
      		   if [ ${decg:0:1} == 0 ]
      		   then decg=`echo $decg | sed -e 's/0//g'`
      		   fi
      	     if [ ${decb:0:1} == 0 ]
      		   then decb=`echo $decb | sed -e 's/0//g'`
      		   fi
      		   let satr=unitr*1000+decr
 	   		   let satg=unitg*1000+decg
      		   let satb=unitb*1000+decb
 		 		   sat=($satr $satg $satg)
 		 	   	 IFS=$'\n'
 		 		   satmax = `echo "${sat[*]}" | sort -nr | head -n1`
 		 		   if (satmax -ge 1000)
 		       then let tv=tv/2
 				   elif (satmax -lt 700)
 			     then let tv=80*tv/satmax
 			     fi
				else echo "Problem with R camera."
						 exit 0
				fi
		 done
		 echo  $tr > Current_R_tint.tmp
		 echo "R integration time: " $tv >> nightmon.log
}


#
# ==================================
#
# main
#
gain=8
max_lum=10000  # 1000000 = 1sec
darkimg="dark-gain16-t100000.dng"
echo  "10000 us" > Current_V_tint.tmp
echo  "10000 us" > Current_R_tint.tmp
# get the site name
/bin/grep "SITE" $HOME/nightmon_config > ligne.tmp
read bidon bidon sitename bidon < ligne.tmp

basepath="/var/www/html/data"
backpath="/home/sand/data"
echo "Johnson V shot"
take_pictureV
y=`date +%Y`
mo=`date +%m`
d=`date +%d`
h=`date +%H`
mi=`date +%M`
s=`date +%S`
basename=`date +%Y-%m-%d_%H-%M-%S`
baseday=`date +%Y-%m-%d`
read  tv toto < Current_V_tint.tmp
# writing to logfile
if [ ! -d $basepath/$y ]
then mkdir $basepath/$y
fi
if [ ! -d $basepath/$y/$mo ]
then /bin/mkdir $basepath/$y/$mo
fi
echo $y $mo $d $h $mi $s " V " $tv $basepath/$y/$m/$basename"_V_"$tv"_"$gain".dng" >> $basepath/$y/$m/nightmon.log
echo "=============================="
# rename pictures
mv capture_1.dng $basepath/$y/$m/$basename"_V_"$tv"_"$gain".dng"
mv capture_1.jpg $basepath/$y/$m/$basename"_V_"$tv"_"$gain".jpg"

echo "Johnson R shot"
take_pictureR
y=`date +%Y`
mo=`date +%m`
d=`date +%d`
h=`date +%H`
mi=`date +%M`
s=`date +%S`
basename=`date +%Y-%m-%d_%H-%M-%S`
read  tr toto < Current_R_tint.tmp
# writing to logfile
if [ ! -d $basepath/$y ]
then mkdir $basepath/$y
fi
if [ ! -d $basepath/$y/$mo ]
then /bin/mkdir $basepath/$y/$mo
fi
echo $y $mo $d $h $mi $s " R " $tr $basename"_R_"$tv"_"$gain".dng" >> $basepath/$y/$m/nightmon.log
echo "=============================="
# rename pictures
mv capture_2.dng $basepath/$y/$m/$basename"_R_"$tv"_"$gain".dng"
mv capture_2.jpg $basepath/$y/$m/$basename"_R_"$tv"_"$gain".jpg"


# check for the night by reading the latest optimal integration time
if [ $tv -lt $max_lum ]
then "Echo too much light. It is probably daytime."
     move_cams.py 2000 1
		 move_cams.py -1500 1
else move_cams.py 2000 1
fi

echo "=============================="
# process sky IMAGES
python3 ProcessNightMon-JVR2H.py -v $basename"_V_"$tv"_"$gain".dng" -r $basename"_R_"$tv"_"$gain".dng" -d $HOME/$darkimg >> $basepath/$y/$m/nightmon.log
# rename pictures
mv Vzeropoint_corr.png $basepath/$y/$m/$basename_Vzeropoint_corr.png
mv Rzeropoint_corr.png $basepath/$y/$m/$basename_Rzeropoint_corr.png
mv VcalSbBkg.png $basepath/$y/$m/$basename_VcalSbBkg.png
mv RcalSbBkg.png $basepath/$y/$m/$basename_RcalSbBkg.png
mv VStars_Match.png $basepath/$y/$m/$basename_VStars_Match.png
mv RStars_Match.png $basepath/$y/$m/$basename_RStars_Match.png
# backup important files
if [ ! -d $backpath/$y ]
then mkdir $backpath/$y
fi
if [ ! -d $backpath/$y/$mo ]
then /bin/mkdir $backpath/$y/$mo
fi
cp -f $basepath/$y/$m/$basename_Vzeropoint_corr.png $backpath/
cp -f $basepath/$y/$m/$basename_Rzeropoint_corr.png $backpath/
cp -f $basepath/$y/$m/$basename_VStars_Match.png $backpath/
cp -f $basepath/$y/$m/$basename_RStars_Match.png $backpath/
cp -f $basepath/$y/$m/nightmon.log $backpath/
cp -f $basepath/$y/$m/"calibrated_"$baseday"_sky.csv"
