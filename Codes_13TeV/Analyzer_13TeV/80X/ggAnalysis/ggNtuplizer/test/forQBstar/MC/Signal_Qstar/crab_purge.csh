#!/bin/tcsh

#foreach coup (f-1p0 f-0p5 f-0p1)
foreach coup (f-0p1 )
#foreach case (500 1000 2000 3000 4000 5000 6000 7000 8000 9000 )
foreach case (500 1000 2000 3000 4000 5000 6000 7000 8000 )

crab purge -d Qstar/crab_job_QstarToGJ_M-${case}_${coup}/ 

end
