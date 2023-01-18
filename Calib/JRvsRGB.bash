B="0 1 2 3 4 5 6 7 8 9 10"
G="0 1 2 3 4 5 6 7 8 9 10"
R=10
for nb in $B
  do for ng in $G
  do cat ./Johnson-RGB.py | sed "s/cecig/$ng/" | sed "s/cecib/$nb/" | sed "s/cecir/$R/" > ./toto.py
    python3 ./toto.py -i 2023-01-11_23-52-07_B_120000000_16.dng  -d /home/aubema/git/nightmon/data/Darks/dark-gain16-t100000.dng -c B -b JR -e 0.05 -m RpiHQ-JFilters  -k stars -z 1 > result.out
    corco=`grep "Correlation coefficient :" result.out | tail -1`
    echo $R $ng $nb $corco >> results.out
  done
done
