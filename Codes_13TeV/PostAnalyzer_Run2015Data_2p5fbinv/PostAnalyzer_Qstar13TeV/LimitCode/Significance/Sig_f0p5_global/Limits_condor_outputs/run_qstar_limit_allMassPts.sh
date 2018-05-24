#!/bin/bash

pwd=$PWD
echo "${pwd}"

File_Dir=/uscms_data/d3/rocky86/13TeV/PostAnalyzer_Qstar13TeV/LimitCode/CMSSW_7_0_9/src/Significance/Sig_f0p5_global

COUNT=6000

for ((i=1;i<=COUNT;i++)); do

coupling=f0p5

###Don't Forget to Complile first

cat>Job_${coupling}_${i}.csh <<EOF
#!/bin/tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh 

cd /uscms_data/d3/rocky86/13TeV/PostAnalyzer_Qstar13TeV/LimitCode/CMSSW_7_0_9/src/
setenv SCRAM_ARCH slc6_amd64_gcc481
cmsenv
#cd ${pwd}
cd \${_CONDOR_SCRATCH_DIR}
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:\${_CONDOR_SCRATCH_DIR}
${File_Dir}/stats > stats_${coupling}_${i}.log
EOF

chmod 775 ${pwd}/Job_${coupling}_${i}.csh

cat>condor_${coupling}_${i} <<EOF
universe = vanilla
Executable = Job_${coupling}_${i}.csh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT

transfer_input_files = ${File_Dir}/stats, ${File_Dir}/Info_ptPhotJet190etaJet2p4dPhi2p5dEta1p8M695LID_76X.root, ${File_Dir}/IP_ptPhotJet190etaJet2p4dPhi2p5dEta1p8M695LID_76X_Qstarf0p5.root, ${BATINSTALLDIR}/lib/libBATmodels.so.3, ${BATINSTALLDIR}/lib/libBAT.so.5, ${BATINSTALLDIR}/lib/libBATmtf.so.0

Queue 1
EOF

condor_submit condor_${coupling}_${i}

done


