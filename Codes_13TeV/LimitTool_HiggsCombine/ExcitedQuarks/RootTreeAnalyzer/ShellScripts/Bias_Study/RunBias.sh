#!/bin/bash

genpdf=dijet  ####dijet,f1,f2,f3,f4,f5
#genpdf=f5  ####dijet,f1,f2,f3,f4,f5
SigStrength=1 ### r=1 corresponds to signal xsec of 10 pb as we have given xsec value of 10 here. 

for fit in f3  #dijet #f5 #f3 f5
do

fitpdf=${fit}  ####dijet,f1,f2,f3,f4,f5

#for mass in 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 
for mass in 1000 
#for mass in 500  ##4500 5000 5500 6000
do

#var0=$(echo "scale=5;${mass}/600" | bc -l)
#var1=$(echo "scale=5;${var0}^-5" | bc -l)
#var2=$(echo "scale=5;(0.22*${var1})+(1*10^-6)" | bc -l)
#bias=${var2}

#var0=$(echo "scale=5;${mass}/500" | bc -l)
#var1=$(echo "e(0.5*l(${var0}))" | bc -l)  ### var0^0.5 = exp(0.5*ln(var0)) doing it another way as bc can not take fractional powers directly.
#var2=$(echo "scale=5;(1.3*${var1})" | bc -l)

#var0=$(echo "scale=30;${mass}^-3" | bc -l)  
#var1=$(echo "scale=30;(7.6*10^10*${var0})" | bc -l)
#bias=$(echo "scale=30;${var1}" | bc -l)

#if (( $(bc <<< "${bias} >= 3.1") )); then
#bias=3.1
#fi

#echo "bias used to decrease pull = ${bias}"

echo "**********************************************"
echo "Running for Mass = ${mass} and fitpdf = ${fit}"
echo "**********************************************"

massPoint=${mass}

xsec=0.01

#outDir=Cards_BiasStudy_Re-miniAOD/${genpdf}_vs_${fitpdf}/BiasResults_${genpdf}_${fitpdf}_r-${SigStrength}_XS-${xsec}_M-${massPoint}/
outDir=Cards_BiasStudy/${genpdf}_vs_${fitpdf}/BiasResults_${genpdf}_${fitpdf}_r-${SigStrength}_M-${massPoint}/
#outDir=Cards_BiasStudy/dijet_vs_${fit}_WithJPSF/BiasResults_${genpdf}_${fitpdf}_r-${SigStrength}_M-${massPoint}/

if [ ! -d ${outDir} ]; then
echo "--------- Making Directory ${outdir} -----------"
mkdir -pv ${outDir}
chmod 775 ${outDir}
fi

##Input fit file will change according to genpdf
#inFitFile=Cards_BiasStudy_Re-miniAOD/${genpdf}BinnedFit/FitResults_ExcitedQuarks2016.root
inFitFile=Cards_BiasStudy/${genpdf}BinnedFit/FitResults_ExcitedQuarks2016.root

intoyFile=${outDir}/mlfitQstar_${massPoint}_lumi-36.813_r-${SigStrength}.000_ExcitedQuarks2016_${genpdf}_${fitpdf}.root
#intoyFile=${outDir}/mlfitQstar_${massPoint}_lumi-35.866_r-${SigStrength}.000_ExcitedQuarks2016_${genpdf}_${fitpdf}.root

python python/RunBias.py -c config/Qstar_bias.config --mass ${massPoint} -m Qstar -d ${outDir} -r ${SigStrength} -l 36.813 --xsec ${xsec} -t 200 --gen-pdf ${genpdf} --fit-pdf ${fitpdf} -b ExcitedQuarks2016 -i ${inFitFile} #--computeBias

python python/PlotBias.py -c config/Qstar_bias.config --mass ${massPoint} -m Qstar -d ${outDir} -r ${SigStrength} -l 36.813 -t ${intoyFile} --gen-pdf ${genpdf} --fit-pdf ${fitpdf} -b ExcitedQuarks2016 #--bias ${bias}

#python python/RunBias.py -c config/Qstar_bias.config --mass ${massPoint} -m Qstar -d ${outDir} -r ${SigStrength} -l 35.866 --xsec ${xsec} -t 200 --gen-pdf ${genpdf} --fit-pdf ${fitpdf} -b ExcitedQuarks2016 -i ${inFitFile}

#python python/PlotBias.py -c config/Qstar_bias.config --mass ${massPoint} -m Qstar -d ${outDir} -r ${SigStrength} -l 35.866 -t ${intoyFile} --gen-pdf ${genpdf} --fit-pdf ${fitpdf} -b ExcitedQuarks2016 #--bias ${bias}

echo "********************************************"
echo "Done for Mass = ${mass} and fitpdf = ${fit}"
echo "********************************************"

done
done














