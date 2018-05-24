#!/bin/tcsh
setenv pwd $PWD

cat>Job_${1}.csh<<EOF
#!/bin/tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh
cd /uscms_data/d3/varun/13TeV/Qstar/CMSSW_7_6_4/src 
setenv SCRAM_ARCH slc6_amd64_gcc493
cmsenv
cd \${_CONDOR_SCRATCH_DIR}
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:\${_CONDOR_SCRATCH_DIR}
./${1}.exe
EOF

chmod 775 ${pwd}/Job_${1}.csh

cat>condor_${1}<<EOF
universe = vanilla
Executable = Job_${1}.csh
Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )
request_disk = 10000000
request_memory = 2100
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = $pwd/${1}.exe,${2}
Output               = ${1}_\$(Cluster)_\$(Process).stdout
Error                = ${1}_\$(Cluster)_\$(Process).stderr
Log                  = ${1}_\$(Cluster)_\$(Process).log
Queue 1
EOF

condor_submit condor_${1}
