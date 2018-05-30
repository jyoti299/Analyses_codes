#!/bin/bash

model=Qstar ##Qstar
coupling=f-1p0
Opti=DEta_Opti  ##DPhi_Opti  ##PhId_Opti  ##JetId_Opti  ##CSVDisc_Opti 

boxName=ExcitedQuarks2016  ##ExcitedQuarks2016  ##Excited1btagQuarks2016

mMin=4600  #1000  #5000 #1000
mMax=5600 #6500 #6000 #5500 #6000

xMin=1e-4 #1e-5
xMax=5e-4 #1e-1 ##5e-4 ##1e-1

python python/PlotOptiLimit.py -m ${model} -b ${boxName} -l 36.813 --massMin ${mMin} --massMax ${mMax} --xsecMin ${xMin} --xsecMax ${xMax} --Opttype ${Opti}  
