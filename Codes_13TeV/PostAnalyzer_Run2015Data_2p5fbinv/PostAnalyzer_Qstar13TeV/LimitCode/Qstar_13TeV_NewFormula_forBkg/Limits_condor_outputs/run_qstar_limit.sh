#!/bin/bash

pwd=$PWD
echo "${pwd}"

File_Dir=/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/LimitCode/CMSSW_7_0_9/src/Qstar_13TeV

#10 mass points
for M in  1000  2000  3000  4000  5000  7000   
#for M in 1000
do 

coupling=f1p0

###Don't Forget to Complile first

cat>Job_${M}_${coupling}.csh <<EOF
#!/bin/tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh 

cd /uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/LimitCode/CMSSW_7_0_9/src/
setenv SCRAM_ARCH slc6_amd64_gcc481
cmsenv
#cd ${pwd}
cd \${_CONDOR_SCRATCH_DIR}
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:\${_CONDOR_SCRATCH_DIR}
${File_Dir}/stats ${M} > stats_${M}_${coupling}.log
EOF

chmod 775 ${pwd}/Job_${M}_${coupling}.csh

cat>condor_${M}_${coupling} <<EOF
universe = vanilla
Executable = Job_${M}_${coupling}.csh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT

transfer_input_files = ${File_Dir}/stats, ${File_Dir}/Histos_data_GJetInvtMass.root, ${File_Dir}/Resonance_Shapes_qstar_f1p0.root, ${BATINSTALLDIR}/lib/libBATmodels.so.3, ${BATINSTALLDIR}/lib/libBAT.so.5, ${BATINSTALLDIR}/lib/libBATmtf.so.0

notify_user = rocky.bala.garg@cern.ch
Output = condor_\$(Cluster)_${M}_${coupling}.stdout
Error = condor_\$(Cluster)_${M}_${coupling}.stderr
Log = condor_\$(Cluster)_${M}_${coupling}.log
Queue 1
EOF

condor_submit condor_${M}_${coupling}

done





