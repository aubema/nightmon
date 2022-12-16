B="-20 -19 -18 -17 -16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
G="-20 -19 -18 -17 -16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
R=20
for nb in $B
  do for ng in $G
  do cat Johnson-RGB.py | sed "s/ceciestrg/$ng/" | sed "s/ceciestrb/$nb/" | sed "s/ceciestrr/$R/" > toto.py
    python3 toto.py -s 2020-09-22_04-22-53_.ARW -d dark-a7s-10s-iso6400.ARW -c A -b JR -m A7S > result.out
    tail -1 result.out >> results.out
  done
done
