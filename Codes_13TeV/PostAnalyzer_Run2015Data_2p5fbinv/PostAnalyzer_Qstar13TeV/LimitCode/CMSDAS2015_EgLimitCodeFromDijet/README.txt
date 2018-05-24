Instructions:
-------------
Please refer to https://twiki.cern.ch/twiki/bin/view/CMS/DijetLimitCode for additional instructions and description.

1) Install BAT (Bayesian Analysis Toolkit):

   wget --no-check-certificate http://www.mppmu.mpg.de/bat/source/BAT-0.9.2.tar.gz
   tar xvzf BAT-0.9.2.tar.gz
   cd BAT-0.9.2
   ./configure --prefix=$HOME
   make
   make install
   cd -

   ### If you are using csh shell
   setenv BATINSTALLDIR    $HOME
   setenv LD_LIBRARY_PATH  {$LD_LIBRARY_PATH}:{$BATINSTALLDIR}/lib
   setenv CPATH            {$BATINSTALLDIR}/include

   ### If you are using Bash shell
   export BATINSTALLDIR=$HOME
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BATINSTALLDIR/lib
   export CPATH=$BATINSTALLDIR/include

   NOTE: BAT needs to be installed only once. However, every time before compiling or using
         the limit code, you will need to initialize the above environment variables.

2) Checkout the limit code package:

   git clone -b CMSDAS2015 git://github.com/CMSDIJET/DijetLimitCode.git LimitCode

   Enter the package directory:

   cd LimitCode

   The package has the following contents:

   Data_and_ResonanceShapes/
   plots/
   Makefile
   README.txt
   binneddata.cc
   binneddata.hh
   fit.cc
   fit.hh
   statistics.cc
   statistics.hh
   stats.cc

   The plots/ subdirectory contains some ROOT plotting macros.

   To start developing your own code, you can start from the stats.cc file.

3) Compile the code:

   make

4) Run the code:

   ./stats MASS

   where MASS is the resonance mass in GeV, for instance

   ./stats 4000
