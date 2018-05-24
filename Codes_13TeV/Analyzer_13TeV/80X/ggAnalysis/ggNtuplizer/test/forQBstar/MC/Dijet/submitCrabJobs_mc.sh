#!/bin/tcsh
#================ Give ALL INPUTS HERE============================
setenv datafile  Datasets_DiJet.txt ##To change
setenv namefile  Outfiles_DiJet.txt ##To change

set events = 20000                                                      # events per job
set tot_events = -1                                                     # total events

#!!set json =   Cert_190456-202305_8TeV_PromptReco_Collisions12_JSON.txt      # Gve the json file to be used
#!!set lumi =       120                   # Lumi per job
#!!set tot_lumi =   -1                   # total LS

setenv user_area Dijet ##To change
setenv user_path /store/user/rocky86/13TeV/Ntuples/80X/MC/${user_area}
setenv out_path /store/user/rgarg/13TeV/Ntuples/80X/MC/${user_area}
setenv cfgDir /uscms_data/d3/rocky86/13TeV/Analyzer_13TeV/80X/CMSSW_8_0_24_patch1/src/ggAnalysis/ggNtuplizer/test
#==================================================================

set file_length1=`wc -l $datafile | cut -c1-2`                  # Count Lines
set file_length2=`wc -l $namefile | cut -c1-2`

set p = 0

#run till the last line of the input files
while ($p != $file_length1)
@ p = ${p} + 1
set Data=`tail -n +$p ${datafile} | head -1`                  # Get name of dataset

#read the name of the output files
set DirName=`tail -n +$p ${namefile} | head -1`              # Get Name of outputfile 

setenv destination_path ${user_path}/${DirName}

if( ! -d ${destination_path} ) then
echo "Making directory ${destination_path}"
eos root://cmseos.fnal.gov mkdir -p ${destination_path}
endif

#set current directory
setenv pwd $PWD

#Give input to the crab.cfg files
cat>crab_${DirName}.py<<EOF
from WMCore.Configuration import Configuration
config = Configuration()
#General
config.section_('General')
config.General.requestName  = 'job_${DirName}'
config.General.workArea     = '${user_area}'
config.General.transferLogs = True
#JobType
config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
#config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles  = ['Spring16_25nsV10_MC.db', 'Spring16_25nsV10_MC_L2Relative_AK8PFchs.txt', 'Spring16_25nsV10_MC_L3Absolute_AK8PFchs.txt']
config.JobType.psetName    = 'run_mc_${DirName}.py'
config.JobType.outputFiles = ['MC_${DirName}.root']
config.JobType.sendExternalFolder = True
#Data
config.section_('Data')
config.Data.inputDataset = '${Data}'
config.Data.splitting    = 'EventAwareLumiBased'
config.Data.unitsPerJob  = ${events}
config.Data.publication  = False
config.Data.totalUnits   = ${tot_events}
#config.Data.ignoreLocality = True
config.Data.outLFNDirBase  ='${out_path}/${DirName}'
#User
config.section_('User')
#Site
config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.blacklist = ["T2_US_Wisconsin"]
#config.Site.whitelist = ['T1_US_FNAL']
#config.Site.whitelist = ['T3_US_FNALLPC']
EOF

#Give the name of outputfile to the PATtuple files--------------------------------------------
sed -e 's|fileName = cms.string('"'"ggtree_mc.root"'"')|fileName = cms.string('"'"MC_${DirName}.root"'"')|' ${cfgDir}/run_mc_80X.py > ${pwd}/run_mc_${DirName}.py

#----------------------------------------------------------------------------------------------
echo "========================================================================================================================"
echo "Submitting Job for ${Data} with Ouput file Name ${DirName}"
echo "========================================================================================================================"

crab submit -c crab_${DirName}.py

if( ! -d junk ) then
echo "Making directory junk"
mkdir junk
chmod 775 -R junk
endif

cp run_mc_${DirName}.py crab_${DirName}.py junk/.
#---------------------------------------------------------------------
end
