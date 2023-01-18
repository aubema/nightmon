B="-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10"
R="-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10"
G=10
for nb in $B
  do for nr in $R
  do cat ./Johnson-RGB.py | sed "s/cecir/$nr/" | sed "s/cecib/$nb/" | sed "s/cecig/$G/" > ./toto.py
    python3 ./toto.py -i 2023-01-11_22-50-07_A_120000000_16.dng  -d /home/aubema/git/nightmon/data/Darks/dark-gain16-t100000.dng -c A -b JV -e 0.102 -m RpiHQ-JFilters  -k stars -z 1 > result.out
    corco=`grep "Correlation coefficient :" result.out | tail -1`
    echo $nr $G $nb $corco >> results.out
  done
done
