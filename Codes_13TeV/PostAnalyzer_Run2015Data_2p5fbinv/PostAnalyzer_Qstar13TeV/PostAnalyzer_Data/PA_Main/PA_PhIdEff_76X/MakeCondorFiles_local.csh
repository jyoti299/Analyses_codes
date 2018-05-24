#!/bin/bash

pwd=$PWD
echo "${pwd}"

cat>Job_${1}.csh <<EOF
#!/bin/tcsh
#source /uscmst1/prod/sw/cms/setup/cshrc prod
source /cvmfs/cms.cern.ch/cmsset_default.csh 
cd /uscms_data/d3/rocky86/13TeV/Analyzer_13TeV/76X/CMSSW_7_6_3_patch2/src/
setenv SCRAM_ARCH slc6_amd64_gcc491
cmsenv
cd \${_CONDOR_SCRATCH_DIR}
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:\${_CONDOR_SCRATCH_DIR} 
${pwd}/${1}.exe
EOF

chmod 775 ${pwd}/Job_${1}.csh

cat>condor_${1} <<EOF
universe = vanilla
Executable = ${pwd}/Job_${1}.csh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = $pwd/${1}.exe
Output = condor_\$(Cluster)_${1}.stdout
Error = condor_\$(Cluster)_${1}.stderr
Log = condor_\$(Cluster)_${1}.log
Queue 1
EOF

condor_submit condor_${1}
