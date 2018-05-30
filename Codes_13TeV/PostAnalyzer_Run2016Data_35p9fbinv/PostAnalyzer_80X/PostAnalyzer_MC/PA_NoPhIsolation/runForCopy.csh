#!/bin/tcsh

setenv InputDir ${1} 
setenv OutputDir ${2}

rm files_to_copy.txt

eosls ${InputDir} > files_to_copy.txt
#ls ${InputDir} > files_to_copy.txt

set count = 0
foreach i (`cat files_to_copy.txt`)

@ count = ${count} + 1
echo "copying ${i}"

#xrdcp root://cmseos.fnal.gov/${InputDir}/${i} ${OutputDir}/${i}
#xrdcp ${InputDir}/${i} root://cmseos.fnal.gov/${OutputDir}/${i} 
xrdcp -r ${InputDir}/${i} root://cmseos.fnal.gov/${OutputDir}/${i} 

end

echo " rTotal ${count} files copied "
rm files_to_copy.txt
 
