#!/bin/bash

PROGRAM=${0##*\/}


MYPROJECT=$1
MYMODEL=$2

function usage() {
cat << EOF
Usage:  $PROGRAM  <project> [<model>]

Creates a new directory named <project> and a new model class named <model>,
i.e., the header and source files for the class, together with a Makefile
for the project.

If <model> is omitted, the name of the model class will be set to <project>.
To only create the project files, use '-' for <model>.
To only create the model class, use '-' for <project>.

EOF
}

if [[ "$1" == '' ]]; then
   usage
   exit 0
fi

if [[ "$2" == "" ]]; then
   MYMODEL=$MYPROJECT
fi

upMYMODEL=`echo $MYMODEL|tr [a-z] [A-Z]`

function m_header() {
cat << EOF
// ***************************************************************
// This file was created using the $PROGRAM script.
// $PROGRAM is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __BAT__|:UP_MODEL:|__H
#define __BAT__|:UP_MODEL:|__H

#include <BAT/BCModel.h>

// This is a |:Model:| header file.
// Model source code is located in file |:Project:|/|:Model:|.cxx

// ---------------------------------------------------------
class |:Model:| : public BCModel
{
   public:

      // Constructors and destructor
      |:Model:|();
      |:Model:|(const char * name);
      ~|:Model:|();

      // Methods to overload, see file |:Model:|.cxx
      void DefineParameters();
      double LogAPrioriProbability(const std::vector<double> &parameters);
      double LogLikelihood(const std::vector<double> &parameters);
      // void MCMCIterationInterface();
};
// ---------------------------------------------------------

#endif

EOF
}

function m_sourcecode() {
cat << EOF
// ***************************************************************
// This file was created using the $PROGRAM script.
// $PROGRAM is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "|:Model:|.h"

#include <BAT/BCMath.h>

// ---------------------------------------------------------
|:Model:|::|:Model:|() : BCModel()
{
   // default constructor
   DefineParameters();
}

// ---------------------------------------------------------
|:Model:|::|:Model:|(const char * name) : BCModel(name)
{
   // constructor
   DefineParameters();
}

// ---------------------------------------------------------
|:Model:|::~|:Model:|()
   // default destructor
{
}

// ---------------------------------------------------------
void |:Model:|::DefineParameters()
{
   // Add parameters to your model here.
   // You can then use them in the methods below by calling the
   // parameters.at(i) or parameters[i], where i is the index
   // of the parameter. The indices increase from 0 according to the
   // order of adding the parameters.

//   AddParameter("x", -10.0, 10.0); // index 0
//   AddParameter("y",  -5.0,  5.0); // index 1
}

// ---------------------------------------------------------
double |:Model:|::LogLikelihood(const std::vector<double> &parameters)
{
   // This methods returns the logarithm of the conditional probability
   // p(data|parameters). This is where you have to define your model.

   double logprob = 0.;

//   double x = parameters.at(0);
//   double y = parameters.at(1);
//   double eps = 0.5;

   // Breit-Wigner distribution of x with nuisance parameter y
//   logprob += BCMath::LogBreitWignerNonRel(x + eps*y, 0.0, 1.0);


   return logprob;
}

// ---------------------------------------------------------
double |:Model:|::LogAPrioriProbability(const std::vector<double> &parameters)
{
   // This method returns the logarithm of the prior probability for the
   // parameters p(parameters).

   double logprob = 0.;

//   double x = parameters.at(0);
//   double y = parameters.at(1);

//   double dx = GetParameter(0)->GetRangeWidth();

//   logprob += log(1./dx);                    // flat prior for x
//   logprob += BCMath::LogGaus(y, 0., 1.0);   // Gaussian prior for y

   return logprob;
}
// ---------------------------------------------------------

EOF
}

function m_makefile1() {
cat << EOF
###################################################################
# This Makefile was created using the $PROGRAM script
# for project $MYPROJECT.
# $PROGRAM is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://www.mppmu.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CFLAGS and LIBS based on the BAT
# installation on your system. Consult the gmake manual for details.
#
###################################################################

# Root variables
ROOTCFLAGS   := \$(shell root-config --cflags)
ROOTLIBS     := -lMinuit \$(shell root-config --libs)

#ROOTLIBS     += -lRooFitCore -lRooFit -lRooStats -lFoam -lMathMore

# compiler and flags
CXX          = g++
CXXFLAGS     =  -g -Wall -fPIC -Wno-deprecated -O2
LD           = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/x86_64-redhat-linux-gnu/bin/ld -m elf_x86_64
LDFLAGS      =  -g -O2
SOFLAGS      = -shared

# standard commands
RM           = rm -f
MV           = mv
ECHO         = echo
CINT         = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/root/5.34.10-cms7//bin/rootcint

# add ROOT flags
CXXFLAGS    += \$(ROOTCFLAGS) 

# ----------------------------------------------------------------------
# The following definitions depend on the setup of the system where
# the project is being compiled. If BAT is installed in the standard
# system search path or the installation directory is defined in the
# BATINSTALLDIR environmental variable then the lines below are correct
# and the compilation will work
CXXFLAGS    += -I. -I./include -I\$(BATINSTALLDIR)/include
LIBS        += -L\$(BATINSTALLDIR)/lib -lBATmodels -lBATmtf -lBAT \$(ROOTLIBS)

#LIBS        +=  -lcuba

# List of all classes (models) used in the program
# Add classes to the end. Backslash indicates continuation
# on the next line
CXXSRCS      = \\
	|:Model:|.cxx

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS      = \$(patsubst %.cxx,%.o,\$(CXXSRCS))
EXEOBJS      =
MYPROGS     = \\
	run|:Project:|

GARBAGE      = \$(CXXOBJS) \$(EXEOBJS) *.o *~ link.d \$(MYPROGS)


# targets
all : project

link.d : \$(patsubst %.cxx,%.h,\$(CXXSRCS))
	\$(CXX) -MM \$(CXXFLAGS) \$(CXXSRCS) > link.d;

-include link.d

%.o : %.cxx
	\$(CXX) \$(CXXFLAGS) -c \$< -o \$@

clean :
	\$(RM) \$(GARBAGE)

project : run|:Project:|.cxx \$(CXXOBJS)
	\$(CXX) \$(CXXFLAGS) -c \$<
	\$(CXX) \$(LDFLAGS) run|:Project:|.o \$(CXXOBJS) \$(LIBS) -o run|:Project:|

print :
	echo compiler  : \$(CXX)
	echo c++ srcs  : \$(CXXSRCS)
	echo c++ objs  : \$(CXXOBJS)
	echo c++ flags : \$(CXXFLAGS)
	echo libs      : \$(LIBS)
	echo so flags  : \$(SOFLAGS)

	echo rootlibs  : \$(ROOTLIBS)

EOF
}

function m_makefile2() {
cat << EOF
###################################################################
# This Makefile was created using the $PROGRAM script
# for project $MYPROJECT.
# $PROGRAM is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://www.mppmu.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CFLAGS and LIBS based on the BAT
# installation on your system. Consult the gmake manual for details.
#
###################################################################

# Root variables
ROOTCFLAGS   := \$(shell root-config --cflags)
ROOTLIBS     := \$(shell root-config --libs) -lMinuit

#ROOTLIBS     += -lRooFitCore -lRooFit -lRooStats -lFoam -lMathMore

# compiler and flags
CXX          = g++
CXXFLAGS     =  -g -Wall -fPIC -Wno-deprecated -O2
LD           = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/x86_64-redhat-linux-gnu/bin/ld -m elf_x86_64
LDFLAGS      =  -g -O2
SOFLAGS      = -shared

# standard commands
RM           = rm -f
MV           = mv
ECHO         = echo
CINT         = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/root/5.34.10-cms7//bin/rootcint

# add ROOT flags
CXXFLAGS    += \$(ROOTCFLAGS) 

# ----------------------------------------------------------------------
# The following definitions depend on the setup of the system where
# the project is being compiled. If BAT is installed in the standard
# system search path or the installation directory is defined in the
# BATINSTALLDIR environmental variable then the lines below are correct
# and the compilation will work
CXXFLAGS    += -I. -I./include -I\$(BATINSTALLDIR)/include
LIBS        += -L\$(BATINSTALLDIR)/lib -lBATmodels -lBATmtf -lBAT \$(ROOTLIBS)

#LIBS        +=  -lcuba

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS      =
EXEOBJS      =
MYPROGS     = \\
	run|:Project:|

GARBAGE      = \$(CXXOBJS) \$(EXEOBJS) *.o *~ link.d \$(MYPROGS)


# targets
all : project

%.o : %.cxx
	\$(CXX) \$(CXXFLAGS) -c \$< -o \$@

clean :
	\$(RM) \$(GARBAGE)

project : run|:Project:|.cxx
	\$(CXX) \$(CXXFLAGS) -c \$<
	\$(CXX) \$(LDFLAGS) run|:Project:|.o \$(LIBS) -o run|:Project:|

print :
	echo compiler  : \$(CXX)
	echo c++ objs  : \$(CXXOBJS)
	echo c++ flags : \$(CXXFLAGS)
	echo libs      : \$(LIBS)
	echo so flags  : \$(SOFLAGS)

	echo rootlibs  : \$(ROOTLIBS)

EOF
}

function m_main1() {
cat << EOF
// ***************************************************************
// This file was created using the $PROGRAM script
// for project $MYPROJECT.
// $PROGRAM is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include "|:Model:|.h"

int main()
{

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

   // create new |:Model:| object
   |:Model:| * m = new |:Model:|();

   BCLog::OutSummary("Test model created");

   // create a new summary tool object
   BCSummaryTool * summary = new BCSummaryTool(m);

   // perform your analysis here

   // normalize the posterior, i.e. integrate posterior
   // over the full parameter space
//   m->Normalize();

   // run MCMC and marginalize posterior wrt. all parameters
   // and all combinations of two parameters
//   m->MarginalizeAll();

   // run mode finding; by default using Minuit
//   m->FindMode();

   // if MCMC was run before (MarginalizeAll()) it is
   // possible to use the mode found by MCMC as
   // starting point of Minuit minimization
//   m->FindMode( m->GetBestFitParameters() );

   // draw all marginalized distributions into a PostScript file
//   m->PrintAllMarginalized("|:Model:|_plots.ps");

   // print all summary plots
//   summary->PrintParameterPlot("|:Model:|_parameters.eps");
//   summary->PrintCorrelationPlot("|:Model:|_correlation.eps");
//   summary->PrintKnowledgeUpdatePlots("|:Model:|_update.ps");

   // calculate p-value
//   m->CalculatePValue( m->GetBestFitParameters() );

   // print results of the analysis into a text file
//   m->PrintResults("|:Model:|_results.txt");

   delete m;
   delete summary;

   BCLog::OutSummary("Test program ran successfully");
   BCLog::OutSummary("Exiting");

   // close log file
   BCLog::CloseLog();

   return 0;

}

EOF
}

function m_main2() {
cat << EOF
// ***************************************************************
// This file was created using the $PROGRAM script
// for project $MYPROJECT.
// $PROGRAM is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

int main()
{

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

   // perform your analysis here

   // close log file
   BCLog::CloseLog();

   BCLog::OutSummary("Test program ran successfully");
   BCLog::OutSummary("Exiting");

   return 0;

}

EOF
}

##########################################################################
##########################################################################
###
###  Main processing is below
###

#
# if we create project
#
if [[ "$MYPROJECT" != "-" ]]; then
   mkdir -p $MYPROJECT

   cat << EOF
--------------------------------------------------------------------------
The new BAT project was created in the directory '$MYPROJECT'.
To test the configuration try to compile the project by running 'make'
inside the directory. In case there are some compilation errors you need
to adjust the parameters inside the 'Makefile'.

Once the program is compiled successfully, you can run it and it should
print some basic information on the screen.

EOF

   if [[ "$MYMODEL" == "-" ]]; then
      m_makefile2 | sed -e "s/|\:Model\:|/$MYMODEL/g; s/|\:Project\:|/$MYPROJECT/g" > $MYPROJECT/Makefile
      m_main2 | sed -e "s/|\:Model\:|/$MYMODEL/g; s/|\:Project\:|/$MYPROJECT/g" > $MYPROJECT/run$MYPROJECT.cxx
   else
      m_header | sed -e "s/|\:UP_MODEL\:|/$upMYMODEL/g; s/|\:Model\:|/$MYMODEL/g; s/|\:Project\:|/$MYPROJECT/g" > $MYPROJECT/$MYMODEL.h
      m_sourcecode | sed -e "s/|\:Model\:|/$MYMODEL/g" > $MYPROJECT/$MYMODEL.cxx
      m_makefile1 | sed -e "s/|\:Model\:|/$MYMODEL/g; s/|\:Project\:|/$MYPROJECT/g" > $MYPROJECT/Makefile
      m_main1 | sed -e "s/|\:Model\:|/$MYMODEL/g; s/|\:Project\:|/$MYPROJECT/g" > $MYPROJECT/run$MYPROJECT.cxx
      echo "Implement your model in files:   $MYMODEL.h  $MYMODEL.cxx"
   fi

cat << EOF
Implement your analysis in file: run$MYPROJECT.cxx

For details consult BAT webpage: http://www.mppmu.mpg.de/bat
--------------------------------------------------------------------------

EOF

#
# if we only create a new model in the current directory
#
else
   m_header | sed -e "s/|\:UP_MODEL\:|/$upMYMODEL/g; s/|\:Model\:|/$MYMODEL/g; s/|\:Project\:|/$MYPROJECT/g" > $MYMODEL.h
   m_sourcecode | sed -e "s/|\:Model\:|/$MYMODEL/g" > $MYMODEL.cxx
   cat << EOF
--------------------------------------------------------------------------
The new BAT model '$MYMODEL' was created in the current directory.

Implement your model in files:   $MYMODEL.h  $MYMODEL.cxx

For details consult BAT webpage: http://www.mppmu.mpg.de/bat
--------------------------------------------------------------------------

EOF

fi

