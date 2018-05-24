#!/bin/tcsh
#source /uscmst1/prod/sw/cms/setup/cshrc prod
source /cvmfs/cms.cern.ch/cmsset_default.csh 
cd /uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/CMSSW_7_4_14/src/
setenv SCRAM_ARCH slc6_amd64_gcc491
cmsenv
cd ${_CONDOR_SCRATCH_DIR}
setenv LD_LIBRARY_PATH /uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/CMSSW_7_4_14/biglib/slc6_amd64_gcc491:/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/CMSSW_7_4_14/lib/slc6_amd64_gcc491:/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/CMSSW_7_4_14/external/slc6_amd64_gcc491/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_14/biglib/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_14/lib/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_14/external/slc6_amd64_gcc491/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/llvm/3.6/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/lib\:${_CONDOR_SCRATCH_DIR} 
/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/PostAnalyzer_Data/PA_Main/PA_JetHT/DATA_JetHT.exe
