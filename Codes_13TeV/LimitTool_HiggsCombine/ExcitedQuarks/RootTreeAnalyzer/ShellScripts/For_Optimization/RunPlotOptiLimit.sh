#!/bin/bash

model=Bstar ##Qstar
coupling=f1p0
Opti=CSVDisc_Opti  ##DEta_Opti  ##DPhi_Opti  ##PhId_Opti  ##JetId_Opti  ##CSVDisc_Opti 

boxName=Excited1btagQuarks2016  ##ExcitedQuarks2016  ##Excited1btagQuarks2016

mMin=1000  #5000 #1000
mMax=5000 #6000 #5500 #6000

xMin=1e-4
xMax=1e-1 ##5e-4 ##1e-1

python python/PlotOptiLimit.py -m ${model} -b ${boxName} -l 24.487 --massMin ${mMin} --massMax ${mMax} --xsecMin ${xMin} --xsecMax ${xMax} --Opttype ${Opti}  
