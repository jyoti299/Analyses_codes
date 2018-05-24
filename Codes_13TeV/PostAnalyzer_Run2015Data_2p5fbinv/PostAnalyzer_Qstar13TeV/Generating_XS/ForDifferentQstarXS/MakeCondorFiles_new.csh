#!/bin/tcsh
setenv pwd $PWD
#echo ${pwd}
cat>Job_${1}_${2}.csh<<EOF
#!/bin/tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh
cd /uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/CMSSW_7_1_14/src/
setenv SCRAM_ARCH slc6_amd64_gcc481
cmsenv
cd \${_CONDOR_SCRATCH_DIR}
cmsRun ${3}
rm XS_M${1}_${2}.root
EOF

chmod 775 ${pwd}/Job_${1}_${2}.csh

cat>condor_${1}_${2}<<EOF
universe = vanilla
Executable = Job_${1}_${2}.csh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = $pwd/${3}
Output               = ${1}_${2}_\$(Cluster)_\$(Process).stdout
Error                = ${1}_${2}_\$(Cluster)_\$(Process).stderr
Log                  = ${1}_${2}_\$(Cluster)_\$(Process).log
Queue 1
EOF

condor_submit condor_${1}_${2}
