
#!/bin/bash

pwd=$PWD
echo "${pwd}"

cat>Job_${1}.csh <<EOF
#!/bin/tcsh
#source /uscmst1/prod/sw/cms/setup/cshrc prod
source /cvmfs/cms.cern.ch/cmsset_default.csh 

cd /uscms_data/d3/rocky86/slc5_amd64_gcc462/Analyzer/CMSSW_5_3_8_patch3/src/
setenv SCRAM_ARCH slc5_amd64_gcc462
cmsenv
cd ${pwd}
${pwd}/${1}.exe
EOF

chmod 775 ${pwd}/Job_${1}.csh

cat>condor_${1} <<EOF
universe = vanilla
Executable = ${pwd}/Job_${1}.csh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Output = condor_\$(Cluster)_${1}.stdout
Error = condor_\$(Cluster)_${1}.stderr
Log = condor_\$(Cluster)_${1}.log
notify_user = rocky.bala.garg@cern.ch
Queue 1
EOF

condor_submit condor_${1}
