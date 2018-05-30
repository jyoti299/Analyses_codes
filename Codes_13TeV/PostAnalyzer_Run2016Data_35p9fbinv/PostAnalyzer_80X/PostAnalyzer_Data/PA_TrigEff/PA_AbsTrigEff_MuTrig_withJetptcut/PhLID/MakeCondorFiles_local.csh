#!/bin/bash

pwd=$PWD
echo "${pwd}"

cat>Job_${1}.csh <<EOF
#!/bin/tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh 
cd /uscms_data/d3/rocky86/13TeV/Analyzer_13TeV/80X/CMSSW_8_0_24_patch1/src
setenv SCRAM_ARCH slc6_amd64_gcc530
cmsenv
cd \${_CONDOR_SCRATCH_DIR}
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:\${_CONDOR_SCRATCH_DIR} 
./${1}.exe
EOF

chmod 775 ${pwd}/Job_${1}.csh

cat>condor_${1} <<EOF
universe = vanilla
Executable = Job_${1}.csh
Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )
request_disk = 1000000
request_memory = 2100
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = $pwd/${1}.exe
Output = condor_\$(Cluster)_${1}.stdout
Error = condor_\$(Cluster)_${1}.stderr
Log = condor_\$(Cluster)_${1}.log
Queue 1
EOF

condor_submit condor_${1}
