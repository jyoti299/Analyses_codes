#!/bin/tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh 
cd /uscms_data/d3/rocky86/13TeV/Analyzer_13TeV/80X/CMSSW_8_0_24_patch1/src
setenv SCRAM_ARCH slc6_amd64_gcc530
cmsenv
cd ${_CONDOR_SCRATCH_DIR}
setenv LD_LIBRARY_PATH /uscms_data/d3/rocky86/13TeV/Analyzer_13TeV/80X/CMSSW_8_0_24_patch1/biglib/slc6_amd64_gcc530:/uscms_data/d3/rocky86/13TeV/Analyzer_13TeV/80X/CMSSW_8_0_24_patch1/lib/slc6_amd64_gcc530:/uscms_data/d3/rocky86/13TeV/Analyzer_13TeV/80X/CMSSW_8_0_24_patch1/external/slc6_amd64_gcc530/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw-patch/CMSSW_8_0_24_patch1/biglib/slc6_amd64_gcc530:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw-patch/CMSSW_8_0_24_patch1/lib/slc6_amd64_gcc530:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw-patch/CMSSW_8_0_24_patch1/external/slc6_amd64_gcc530/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_24/biglib/slc6_amd64_gcc530:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_24/lib/slc6_amd64_gcc530:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/llvm/3.8.0-giojec2/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/lib\:${_CONDOR_SCRATCH_DIR} 
./singlePhotonRun2016Bv2_401_450.exe
