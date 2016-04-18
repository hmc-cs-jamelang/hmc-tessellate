#!/bin/sh
rm testdata.txt
touch testdata.txt
rm testraw.txt
touch testraw.txt

for i in 10 20 40 70 100 200 400 700 1000 2000 4000 7000 10000 20000 40000 70000 100000 200000 400000 700000 1000000 2000000 4000000 7000000 10000000 # 20000000 40000000 70000000 100000000
do
	make test n=5000 2>/dev/null 1>&2 && make test n=$i 2>> testraw.txt
done

sed -n '/ms/p' testraw.txt > temp
echo "voro++ data:" >> testdata.txt
sed -n '1~2p' temp >> testdata.txt
echo "hmc-tessellate data:" >> testdata.txt
sed -n '0~2p' temp >> testdata.txt
rm temp
