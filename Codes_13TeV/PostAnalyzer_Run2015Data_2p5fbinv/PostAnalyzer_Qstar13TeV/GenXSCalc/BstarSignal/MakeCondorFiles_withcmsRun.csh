#!/bin/bash

pwd=$PWD
echo "${pwd}"

ana_Dir=/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/GenXSCalc/CMSSW_7_5_0/src

cat>Job.csh <<EOF
#!/bin/tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh 
cd /uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/GenXSCalc/CMSSW_7_5_0/src/
setenv SCRAM_ARCH slc6_amd64_gcc491
cmsenv
cd \${_CONDOR_SCRATCH_DIR}
${pwd}/run_ana.sh ${1}
EOF

chmod 775 ${pwd}/Job.csh

cat>condor.jdl <<EOF
universe = vanilla
Executable = ${pwd}/Job.csh
Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )
request_disk = 10000000
request_memory = 2100
use_x509userproxy = true
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = ${ana_Dir}/das_client.py, ${ana_Dir}/ana.py, ${pwd}/run_ana.sh, ${pwd}/${1}
Output = condor_\$(Cluster).stdout
Error = condor_\$(Cluster).stderr
Log = condor_\$(Cluster).log
Queue 1
EOF

condor_submit condor.jdl
