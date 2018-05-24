#!/bin/bash

pwd=$PWD
echo "${pwd}"

File_Dir=/uscms_data/d3/rocky86/13TeV/PostAnalyzer_Qstar13TeV/LimitCode/CMSSW_7_0_9/src/Significance/Sig_f0p1_global

COUNT=500

for ((i=1;i<=COUNT;i++)); do

Inseed=7629
JobID=$[Inseed*i]

for M in  1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000 4100 4200 4300 4400 4500 4600 4700 4800 4900 5000
do 

coupling=f0p1

###Don't Forget to Complile first

cat>Job_${M}_${coupling}_${i}.csh <<EOF
#!/bin/tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh 

cd /uscms_data/d3/rocky86/13TeV/PostAnalyzer_Qstar13TeV/LimitCode/CMSSW_7_0_9/src/
setenv SCRAM_ARCH slc6_amd64_gcc481
cmsenv
#cd ${pwd}
cd \${_CONDOR_SCRATCH_DIR}
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:\${_CONDOR_SCRATCH_DIR}
${File_Dir}/stats ${M} qq ${JobID} > stats_${M}_${coupling}_${i}.log
EOF

chmod 775 ${pwd}/Job_${M}_${coupling}_${i}.csh

cat>condor_${M}_${coupling}_${i} <<EOF
universe = vanilla
Executable = Job_${M}_${coupling}_${i}.csh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT

transfer_input_files = ${File_Dir}/stats, ${File_Dir}/Info_ptPhotJet190etaJet2p4dPhi2p5dEta1p8M695LID_76X.root, ${File_Dir}/IP_ptPhotJet190etaJet2p4dPhi2p5dEta1p8M695LID_76X_Qstarf0p1.root, ${BATINSTALLDIR}/lib/libBATmodels.so.3, ${BATINSTALLDIR}/lib/libBAT.so.5, ${BATINSTALLDIR}/lib/libBATmtf.so.0

Queue 1
EOF

condor_submit condor_${M}_${coupling}_${i}

done



done




