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
		 echo "Taking pictures"
		 read tv bidon < Current_V_tint.tmp
		 read tr bidon < Current_R_tint.tmp
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
		 while [ $satmax -gt 75 ]
		 do	rm -f capture_1*
		 		python3 captureA.py -i $tv
				if [ test -f capture_1.dng ]
		    then lisc perc capture_1.dng -p 99.9 | sed -e 's/=//g' | sed -e 's/R//g' | sed -e 's/G//g' | sed -e 's/B//g'| sed -e 's/\./ /g' > capture.tmp
		 	   	 read unitr decr unitg decg unitb decb bidon  < capture.tmp
     		   # remove leading zero to the sky brightness
     		   if [ ${unitr:0:1} == 0 ]
     		   then unitr=`echo $unitr | sed -e 's/0//g'`
     		   fi
     		   if [ ${unitg:0:1} == 0 ]
     		   then unitg=`echo $unitg | sed -e 's/0//g'`
     		   fi
     	     if [ ${unitb:0:1} == 0 ]
     		   then unitb=`echo $unitb | sed -e 's/0//g'`
     		   fi
     		   let satr=unitr*100+decr
	   		   let satg=unitg*100+decg
     		   let satb=unitb*100+decb
		 		   sat=($satr $satg $satg)
		 	   	 IFS=$'\n'
		 		   satmax = `echo "${sat[*]}" | sort -nr | head -n1`
		 		   if (satmax -ge 100)
		       then let tv=tv/2
		 		   elif (satmax -lt 75)
			     then let tv=75*tv/satmax
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
		 while [ $satmax -gt 75 ]
		 do	rm -f capture_2*
		 		python3 captureB.py -i $tr
				if [ test -f capture_2.dng ]
				then lisc perc capture_2.dng -p 99.9 | sed -e 's/=//g' | sed -e 's/R//g' | sed -e 's/G//g' | sed -e 's/B//g'| sed -e 's/\./ /g' > capture.tmp
				 		read unitr decr unitg decg unitb decb bidon  < capture.tmp
		     		# remove leading zero to the sky brightness
		     		if [ ${unitr:0:1} == 0 ]
		     		then unitr=`echo $unitr | sed -e 's/0//g'`
		     		fi
		     		if [ ${unitg:0:1} == 0 ]
		     		then unitg=`echo $unitg | sed -e 's/0//g'`
		     		fi
		     		if [ ${unitb:0:1} == 0 ]
		     		then unitb=`echo $unitb | sed -e 's/0//g'`
		     		fi
		     		let satr=unitr*100+decr
			   		let satg=unitg*100+decg
		     		let satb=unitb*100+decb
				 		sat=($satr $satg $satg)
				 		IFS=$'\n'
				 		satmax = `echo "${sat[*]}" | sort -nr | head -n1`
				 		if (satmax -ge 100)
				    	then let tr=tr/2
				 		elif (satmax -lt 75)
					  	then let tr=75*tr/satmax
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
max_lum=10000  # 1000000 = 1sec
# get the site name
/bin/grep "SITE" $HOME/nightmon_config > ligne.tmp
read bidon bidon sitename bidon < ligne.tmp
# Sets gpio 13 as an output for the LED
#=====
#
# main loop
#
time1=`date +%s`
i=0
y=`date +%Y`
mo=`date +%m`
d=`date +%d`
h=`date +%H`
mi=`date +%M`
s=`date +%S`
basename=`date +%Y-%m-%d_%H-%M-%S`

#basepath="/var/www/html/data"
basepath="./test"


take_pictureV
y=`date +%Y`
mo=`date +%m`
d=`date +%d`
h=`date +%H`
mi=`date +%M`
s=`date +%S`
basename=`date +%Y-%m-%d_%H-%M-%S`
read  tv toto < Current_V_tint.tmp
# writing to logfile
if [ ! -d $basepath/$y ]
then mkdir $basepath/$y
fi
if [ ! -d $basepath/$y/$mo ]
then /bin/mkdir $basepath/$y/$mo
fi
echo $y $mo $d $h $mi $s " V " $tv $basepath/$y/$m/V-$tv-$basename.dng >> $basepath/$y/$m/nightmon.log
echo "=============================="
# rename pictures
cp capture_1.dng $basepath/$y/$m/$basename_V.dng
cp capture_1.jpg $basepath/$y/$m/$basename_V.jpg


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
echo $y $mo $d $h $mi $s " R " $tr $basepath/$y/$m/R-$tr-$basename.dng >> $basepath/$y/$m/nightmon.log
echo "=============================="
# rename pictures
cp capture_1.dng $basepath/$y/$m/$basename_V.dng
cp capture_1.jpg $basepath/$y/$m/$basename_V.jpg


# check for the night by reading the latest optimal integration time
if [ $tv -lt $max_lum ]
then "Echo too much light. It is probably daytime."
     move_cams.py 2000 1
		 move_cams.py -1500 1
else move_cams.py 2000 1
fi

echo "=============================="
# process sky IMAGES
python3 ProcessNightMon-JVR2H.py -v capture_1.dng -r capture_2.dng >> $basepath/$y/$m/nightmon.log
# rename pictures
mv BackgroundV.npy $basepath/$y/$m/$basename_BackgroundV.npy
mv BackgroundR.npy $basepath/$y/$m/$basename_BackgroundR.npy
mv Vzeropoint_corr.png $basepath/$y/$m/$basename_Vzeropoint_corr.png
mv Rzeropoint_corr.png $basepath/$y/$m/$basename_Rzeropoint_corr.png
mv VcalSbBkg.png $basepath/$y/$m/$basename_VcalSbBkg.png
mv RcalSbBkg.png $basepath/$y/$m/$basename_RcalSbBkg.png
mv VStars_Match.png $basepath/$y/$m/$basename_VStars_Match.png
mv RStars_Match.png $basepath/$y/$m/$basename_RStars_Match.png
