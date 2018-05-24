#!/bin/tcsh

setenv coupling f0p1
setenv coup f-0p1
foreach case (500 1000 1500 2000 2500 3000 3500 4000 4500 5000 )

crab report -d ${coupling}/Signal_Bstar/${coupling}/crab_job_BstarToGJ_M-${case}_${coup}/ 

end
