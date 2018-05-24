-------------
Instructions:
-------------
Please refer to https://twiki.cern.ch/twiki/bin/view/CMS/DijetLimitCode for additional instructions and description.


1) Set up your CMSSW working area:

   setenv SCRAM_ARCH slc6_amd64_gcc481
   cmsrel CMSSW_7_0_9
   cd CMSSW_7_0_9/test
   cmsenv

   NOTE: You can skip this step if you already have your working area set up. Nevertheless,
         you will still need to initialize the CMSSW environment by calling 'cmsenv'.

2) Install BAT (Bayesian Analysis Toolkit):

   wget --no-check-certificate http://www.mppmu.mpg.de/bat/source/BAT-0.9.2.tar.gz
   tar xvzf BAT-0.9.2.tar.gz
   cd BAT-0.9.2
   ./configure --prefix=$HOME
   make
   make install
   cd -

   setenv BATINSTALLDIR    $HOME
   setenv LD_LIBRARY_PATH  {$LD_LIBRARY_PATH}:{$BATINSTALLDIR}/lib
   setenv CPATH            {$BATINSTALLDIR}/include

   Alternatively, instead of executing the above three lines, you can just source BAT_init.csh

   source ./BAT_init.csh

   NOTE: BAT needs to be installed only once. However, every time before compiling or using
         the limit code, you will need to initialize the above environment variables.

3) Checkout the limit code package:

   git clone git://github.com/CMSDIJET/DijetLimitCode.git DijetLimitCode

   Enter the package directory:

   cd LimitCode

   The package has the following contents:

   Data_and_ResonanceShapes/
   plots/
   jobs/
   BAT_init.csh
   Makefile
   README.txt
   binneddata.cc
   binneddata.hh
   fit.cc
   fit.hh
   statistics.cc
   statistics.hh
   stats.cc

   The plots/ subdirectory contains some ROOT plotting macros and scripts and the jobs/ subdirectory
   contains some Condor job submission scripts.

   To start developing your own code, you can start from the stats.cc file.

4) Compile the code:

   make

5) Run the code:

   ./stats MASS FINAL_STATE

   where MASS is the resonance mass in GeV and FINAL_STATE is the type of two-parton final state, for instance

   ./stats 4000 qg

--------------------
Running using Condor:
--------------------

In order to speed up the process of running the limit code for multiple mass points, it is possible to use
the Condor batch queue. A Python script is provided that can automatically create and submit jobs for multiple
mass points and different final states.

To create and submit Condor jobs, run:

jobs/submit_condor_jobs.py FINAL_STATE MASS_MIN MASS_MAX STEP NPES NJOBS

where FINAL_STATE is the type of two-parton final state, MASS_MIN is the initial mass point, MASS_MAX is the final
mass point, and STEP defines the step size between the initial and final mass points. All mass-related quantities
are defined in units of GeV. The final two input arguments, NPES and NJOBS, are optional. The NPES argument
defines the number of pseudo-experiments per job and NJOBS defines the number of jobs per mass point. The total
number of pseudo-experiments per mass point will be NPES*NJOBS. If NPES is not specified, the hardcoded value
will be used. If NJOBS is not specified, only one job per mass point will be run. A possible example is

jobs/submit_condor_jobs.py qg 1500 7200 100 20 10

which will run limits for qg resonances starting from a mass point at 1.2 TeV going up to 7 TeV in steps of 100 GeV.
10 jobs per mass point will be created each running 20 pseudo-experiments resulting in a total of 200 pseudo-experiments
per mass point.

To combine the output from multiple jobs and get only one .log and .root file per mass point, run:

jobs/combine_output.py FINAL_STATE MASS_MIN MASS_MAX STEP

or for the above example

jobs/combine_output.py qg 1500 7200 100

Finally, to clean up your working directly from Condor files that are no longer needed, execute:

rm run_limits_*.jdl condor_*.log condor_*.stdout condor_*.stderr
