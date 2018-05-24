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

