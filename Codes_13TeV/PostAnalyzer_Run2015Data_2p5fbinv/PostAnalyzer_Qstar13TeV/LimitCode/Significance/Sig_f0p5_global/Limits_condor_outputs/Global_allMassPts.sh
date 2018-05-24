#!/bin/bash


#Making Dir and moving all the stats.log files into this
FILES_DIR=STATS_LOG

if [ ! -d ${FILES_DIR} ]; then
echo "Making Directory ${FILES_DIR}"
mkdir ${FILES_DIR}
chmod 775 ${FILES_DIR}
fi

echo "MOVING ALL STATS LOG FILES TO STATS_LOG DIR"
mv stats*.log ${FILES_DIR}/
#-------------------------------------------------------

#Finding all the stats.log files with Fit status: FAILED and deleting them
for f in STATS_LOG/*
do

grep "Fit status: FAILED" "$f" > L

#-s gives size of file
if [ -s L ]; then
echo "Deleting $f: FAILED FIT STATUS" 
rm "$f"
fi

rm L

done
#--------------------------------------------------------------------------


#Making Dirs
Dir1=STATS_OUT

if [ ! -d ${Dir1} ]; then
echo "Making Directory ${Dir1}"
mkdir ${Dir1}
chmod 775 ${Dir1}
fi

Dir2=STATS_SORT

if [ ! -d ${Dir2} ]; then
echo "Making Directory ${Dir2}"
mkdir ${Dir2}
chmod 775 ${Dir2}
fi
#-----------

for file in STATS_LOG/*
do

filename=$(echo "$file" | awk -F'[_.]' '{print $4}')

echo "Getting Significance: ${filename}"

#Moving all the Significance numbers in different files
grep "Significance(data)" "$file" > ${Dir1}/${filename}

#sorting all the significance numbers 
echo "Sorting Significance: ${filename}"
sort -nrk6 ${Dir1}/${filename} > ${Dir2}/${filename}

#Taking out the max significance numbers in a single file
echo "Getting Maximum Significance: ${filename}"
head -n 1 ${Dir2}/${filename} >> Max_Sig.txt

done

#finally sorting to get the max numbers on the top
echo "Final Sorting"
sort -nrk6 Max_Sig.txt > Max_Sig_sorted.txt

