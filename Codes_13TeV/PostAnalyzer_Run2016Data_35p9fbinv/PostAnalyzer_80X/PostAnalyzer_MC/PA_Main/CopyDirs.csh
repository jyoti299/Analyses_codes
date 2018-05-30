#!/bin/tcsh

setenv InputDir ${1} 
setenv OutputDir ${2}

eosls ${InputDir} > dirs_to_copy.txt

setenv datafile dirs_to_copy.txt
set file_length=`wc -l ${datafile} | cut -c1-2` 

set p = 0
#run till the last line of the input files      
while ($p != $file_length)
@ p = ${p} + 1
set Dir=`tail -n +$p ${datafile} | head -1`                  # Get name of dataset                  

setenv destination_path ${OutputDir}/${Dir}

if( ! -d ${destination_path} ) then
echo "Making directory ${destination_path}"
eos root://cmseos.fnal.gov mkdir -p ${destination_path}
endif

eosls ${InputDir}/${Dir} > files_to_copy.txt

set count = 0
foreach i (`cat files_to_copy.txt`)

@ count = ${count} + 1
echo "copying ${i}"

#xrdcp root://cmseos.fnal.gov/${InputDir}/${i} ${OutputDir}/${i}
#xrdcp ${InputDir}/${i} root://cmseos.fnal.gov/${OutputDir}/${i} 
xrdcp  root://cmseos.fnal.gov/${InputDir}/${Dir}/${i} root://cmseos.fnal.gov/${OutputDir}/${Dir}/${i} 

end

echo " rTotal ${count} files copied "
rm files_to_copy.txt
 
ehco "Done for ${Dir}"

end
