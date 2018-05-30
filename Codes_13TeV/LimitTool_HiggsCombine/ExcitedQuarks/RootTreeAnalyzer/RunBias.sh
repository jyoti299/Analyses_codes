#!/bin/bash

model=Qstar
boxName=ExcitedQuarks2016

#genpdf=dijet  ####dijet,f1,f2,f3,f4,f5
fitpdf=dijet  ####dijet,f1,f2,f3,f4,f5
SigStrength=1 ### r=1 corresponds to signal xsec of 10 pb as we have given xsec value of 10 here. 

#for fit in f1 f2 #f5 #f3 f5
for gen in f6 #f5 #f3 f5
do

genpdf=${gen}  ####dijet,f1,f2,f3,f4,f5

#for mass in 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 
for mass in 2000  
do

nToys=100

echo "**********************************************"
echo "Running for Mass = ${mass} and genpdf = ${genpdf}"
echo "**********************************************"

massPoint=${mass}

xsec=0.004

#outDir=Cards_BiasStudy/${genpdf}_vs_${fitpdf}/BiasResults_${genpdf}_${fitpdf}_r-${SigStrength}_XS-${xsec}_M-${massPoint}/
outDir=Cards_BiasStudy_Inverted/${model}/${genpdf}_vs_${fitpdf}/BiasResults_${genpdf}_${fitpdf}_r-${SigStrength}_XS-${xsec}_M-${massPoint}/

if [ ! -d ${outDir} ]; then
echo "--------- Making Directory ${outdir} -----------"
mkdir -pv ${outDir}
chmod 775 ${outDir}
fi

config=config/Qstar_bias.config

##Input fit file will change according to genpdf
inFitFile=Cards_BiasStudy_Inverted/${model}/${genpdf}BinnedFit/FitResults_${boxName}.root
#inFitFile=Cards_GOF/${model}/F6BinnedFit/FitResults_${boxName}.root

#intoyFile=${outDir}/mlfit${model}_${massPoint}_lumi-35.866_r-${SigStrength}.000_${boxName}_${genpdf}_${fitpdf}.root
intoyFile=${outDir}/mlfit${model}_${massPoint}_lumi-36.813_r-${SigStrength}.000_${boxName}_${genpdf}_${fitpdf}.root

python python/RunBias.py -c ${config} --mass ${massPoint} -m ${model} -d ${outDir} -r ${SigStrength} -l 36.813 --xsec ${xsec} -t ${nToys} --gen-pdf ${genpdf} --fit-pdf ${fitpdf} -b ${boxName} -i ${inFitFile}
#python python/RunBias.py -c ${config} --mass ${massPoint} -m ${model} -d ${outDir} -r ${SigStrength} -l 35.866 --xsec ${xsec} -t ${nToys} --gen-pdf ${genpdf} --fit-pdf ${fitpdf} -b ${boxName} -i ${inFitFile}

python python/PlotBias.py -c ${config} --mass ${massPoint} -m ${model} -d ${outDir} -r ${SigStrength} -l 36.813 -t ${intoyFile} --gen-pdf ${genpdf} --fit-pdf ${fitpdf} -b ${boxName} #--bias ${bias}
#python python/PlotBias.py -c ${config} --mass ${massPoint} -m ${model} -d ${outDir} -r ${SigStrength} -l 35.866 -t ${intoyFile} --gen-pdf ${genpdf} --fit-pdf ${fitpdf} -b ${boxName} #--bias ${bias}

echo "********************************************"
echo "Done for Mass = ${mass} and genpdf = ${genpdf}"
echo "********************************************"

done
done














