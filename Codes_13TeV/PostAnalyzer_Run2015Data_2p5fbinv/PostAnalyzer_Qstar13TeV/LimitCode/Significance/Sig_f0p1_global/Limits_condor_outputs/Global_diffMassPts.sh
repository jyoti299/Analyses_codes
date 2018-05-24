#!/bin/bash

f=f0p1
FILES_DIR=STATS_LOG

if [ ! -d ${FILES_DIR} ]; then
echo "Making Directory ${FILES_DIR}"
mkdir ${FILES_DIR}
chmod 775 ${FILES_DIR}
fi

echo "MOVING ALL STATS LOG FILES TO STATS_LOG DIR"
mv stats*.log ${FILES_DIR}/

Final_Dir=STATS_FINAL

if [ ! -d ${Final_Dir} ]; then
echo "Making Directory ${Final_Dir}"
mkdir ${Final_Dir}
chmod 775 ${Final_Dir}
fi

for M in 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500
do

if [ ! -d ${M} ]; then
echo "Making Directory ${M}"
mkdir ${M}
chmod 775 ${M}
fi

cp ${FILES_DIR}/stats_${M}_${f}_*.log ${M}/
grep "Significance" ${M}/stats_${M}_${f}_*.log >> ${M}/stats_${M}.log
awk '{ print $3 }'  ${M}/stats_${M}.log >> ${M}/stats_${M}_awk.log
sort -n -r ${M}/stats_${M}_awk.log >> ${M}/stats_${M}_sorted.log

cp ${M}/stats_${M}_sorted.log ${Final_Dir}/

rm -r ${M}

done
