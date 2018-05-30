------------
Instructions
------------
1) set up CMSSW working area:

setenv SCRAM_ARCH slc6_amd64_gcc491
cmsrel CMSSW_7_4_14
cd CMSSW_7_4_14/src
cmsenv

2) Download latest HiggsCombine version:

git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v6.3.1
scramv1 b clean; scramv1 b # always make a clean build.


3) Checking out the package:


4) Details:
Limit scripts have been made by taking the idea from Dijet limit code present in CMSDIJET/DijetRootTreeAnalyzer/ and DijetShapeInterpolator/

Main scipts, to use, are present at ExcitedQuarks/RootTreeAnalyzer/, ExcitedQuarksShapeInterpolator/ and Fisher_test/ 

The HiggsAnalysis_modified/CombinedLimit/ area is the HiggsCombine package but with some added files in src/ and interface/ areas for additional background functions.
The procedure to add a polynomial function, for it to be used by the limit tool, is to make a .cc file in src/ directory and a .h file in interface/ directory with the
information of the function (see for example src/RooDijetBinPdf.cc and interface/RooDijetBinPdf.h) and add the information of the function in src/classes.h
and src/classes_def.xml files. 

If you want to add an additional function to use, then take the idea from HiggsAnalysis_modified/CombinedLimit/ and make yours. If you want to use one already present
in HiggsAnalysis_modified/CombinedLimit/ area, just copy .cc and .h files in respective areas and add information in src/classes.h and src/classes_def.xml files.
You can copy the whole HiggsAnalysis_modified/CombinedLimit/ area to HiggsAnalysis/CombinedLimit/ and use that as well if the latest version is not very different from
v6.3.1.  



