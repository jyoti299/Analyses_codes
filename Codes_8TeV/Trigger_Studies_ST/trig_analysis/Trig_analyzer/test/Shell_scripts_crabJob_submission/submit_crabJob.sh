#!/bin/bash



outFile=trig_crab_outputFile.txt
tagFile=trig_tag.txt

outFile_length=$(wc -l $(outFile) | cut -c1-2)
tagFile_length=$(wc -l $(tagFile) | cut -c1-2)

if[ $outFile_length -eq $tagFile_length ]; then
echo "No. of lines in both outfile and tagfile are $outFile_length"

p=0

while[ $p -ne $outFile_length ]; do
((p++))

outputFile=$(tail -n +$p $(outFile) | head -1)

tagName=$(tail -n +$p $(tagFile) | head -1)




cat<crab_$(outputFile).cfg <<EOF


[CRAB]
jobtype    = cmssw
scheduler  = condor
#scheduler  = remoteGlidein
#use_server = 0


[CMSSW]
datasetpath             = /DoubleElectron/Run2012D-ZElectron-22Jan2013-v1/RAW-RECO
#output_file             = Standard_HLTLevel_corrections.root
output_file            = Laser_corrections_perCrystal_1day.root
#pset                    = trig_analyzer_cfg.py
pset                   = hltconfig.py
lumi_mask               = Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt
total_number_of_lumis   = -1
lumis_per_job           = 10
#events_per_job         = 20000
#total_number_of_events = -1
#number_of_jobs         =
get_edm_output          = 1

[GRID]
#dont_check_proxy = 1
#rb                     = CERN
#proxy_server           = myproxy.cern.ch
virtual_organization   = cms
#retry_count            = 0
#ce_black_list=dc2-grid-70.brunel.ac.uk, grid01.physics.uoi.gr, grcreamce01.inr.troitsk.ru, T2_CN_Beijing, T2_CN_Beijing, T2_RU_RRC_KI, oberon.hep.kbf\
i.ee, mars.hep.kbfi.ee, cream1.hep.kbfi.ee$
#se_white_list         =

[USER]
#eMail                  = rocky.bala.garg@cern.ch

###To store the output locally
return_data           = 0
#outputdir = /uscms_data/d3/rocky86/slc5_amd64_gcc462/Trigger_studies/CMSSW_5_3_8_patch3/src/trig_analysis/Trig_analyzer/test/HLT_laser_corrections/St\
andard_corrections_

###To store output in eos or dcache
copy_data             = 1
###storage element and storage path for eos
storage_element      = cmseos.fnal.gov
storage_path         = /srm/v2/server?SFN=/eos/uscms/store/user/rocky86
###storage element and storage path for dcache
#storage_element       = cmssrm.fnal.gov
#storage_path          = /srm/managerv2?SFN=/11/store/user/rocky86
user_remote_dir       = /HLT_laser_corrections/PerCrystal_1day_corrections
ui_working_dir        = Crab_perCrystal_1day_corrections

EOF



done

fi
