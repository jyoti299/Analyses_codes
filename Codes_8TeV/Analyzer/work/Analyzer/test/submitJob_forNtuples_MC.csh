#!/bin/tcsh
#================ Give ALL INPUTS HERE============================
setenv datafile  gjdataset.txt
setenv namefile  gjoutputfile.txt
set events = 5000                                                         # events per job
set tot_events = -1                                                        # total events
set json =   Cert_190456-202305_8TeV_PromptReco_Collisions12_JSON.txt      # Gve the json file to be used
set lumi =       120                   # Lumi per job
set tot_lumi =   -1                   # total LS
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


#setenv destination_path /eos/uscms/store/user/varun/Summer12/53X/PhotonJet/
setenv destination_path  /pnfs/cms/WAX/11/store/user/varun/2012/MC/Summer12/53X/QstarSignal/

setenv destination_dir ${destination_path}${DirName}

if( ! -d ${destination_path} ) then
echo "Making directory ${destination_path}"
mkdir ${destination_path}
chmod 775 -R ${destination_path}
endif

if( ! -d ${destination_dir} ) then
echo "Making directory ${destination_dir}"
mkdir ${destination_path}${DirName}
chmod 775 -R ${destination_path}${DirName}
endif
#set direcotry name in store
setenv filedir /2012/MC/Summer12/53X/QstarSignal/${DirName}

#set current directory
setenv pwd $PWD

#Give input to the crab.cfg files
cat>crab_${DirName}.cfg<<EOF
[CMSSW]
# FOR DATA to select good runs only
#lumi_mask             = ${json}      
#total_number_of_lumis = ${tot_lumi} 
#lumis_per_job         = ${lumi}    
events_per_job         = ${events}
total_number_of_events = ${tot_events}
#number_of_jobs         =
pset                   = config_53X_MC_Sig_cfg.py 
datasetpath            = ${Data}
#runselection          =
output_file            = AOD_Output_${DirName}.root

[USER]
eMail                  = varun.sharma@cern.ch
return_data            = 0
copy_data              = 1
ui_working_dir         = ${DirName}

#storage_element        = cmseos.fnal.gov
storage_element        = cmssrm.fnal.gov

## and the SE directory (or the mountpoint) that has to be writable from all
#### LNL SRM
#storage_path           = /srm/v2/server?SFN=/eos/uscms/store/user/varun
storage_path           = /srm/managerv2?SFN=/11/store/user/varun
user_remote_dir        = ${filedir}
publish_data           = 0

[GRID]
#rb                     = CERN
#proxy_server           = myproxy.cern.ch
virtual_organization   = cms
retry_count            = 0
#ce_black_list=dc2-grid-70.brunel.ac.uk, grid01.physics.uoi.gr, grcreamce01.inr.troitsk.ru, T2_CN_Beijing, T2_RU_RRC_KI, oberon.hep.kbfi.ee, mars.hep.kbfi.ee, cream1.hep.kbfi.ee 
#se_white_list         =

[CRAB]
scheduler              = condor
#scheduler              = remoteGlidein
jobtype                = cmssw
#use_server             = 1
EOF

#Give the name of outputfile to the PATtuple files--------------------------------------------
sed -e 's|outFile          = cms.untracked.string|outFile          = cms.untracked.string('"'"AOD_Output_${DirName}.root"'"'),|'   ${pwd}/config_53X_MC_Sig_cfg.py > ${pwd}/config_53X_MC_Sig_cfg_tmp.py

cp config_53X_MC_Sig_cfg.py config_53X_MC_Sig_cfg_orig.py
mv config_53X_MC_Sig_cfg_tmp.py config_53X_MC_Sig_cfg.py
#----------------------------------------------------------------------------------------------
echo "========================================================================================================================"
echo "Submitting Job for ${Data} with Ouput file Name ${DirName}"
echo "========================================================================================================================"

#sleep for 10 minuts

#create and submit the job---------------------------------------------
if(-f log_Submit_${DirName})then
crab -create -submit -cfg crab_${DirName}.cfg >> log_Submit_${DirName}
else
crab -create -submit -cfg crab_${DirName}.cfg > log_Submit_${DirName}
endif
#--------------------------------------------------------------------


if( ! -d junk ) then
echo "Making directory junk"
mkdir junk
chmod 775 -R junk
endif


#now put back the original file at its place
cp config_53X_MC_Sig_cfg.py junk/config_53X_MC_Sig_${DirName}_cfg.py
mv config_53X_MC_Sig_cfg_orig.py config_53X_MC_Sig_cfg.py 
cp crab_${DirName}.cfg junk/crab_${DirName}.cfg
#---------------------------------------------------------------------

end
