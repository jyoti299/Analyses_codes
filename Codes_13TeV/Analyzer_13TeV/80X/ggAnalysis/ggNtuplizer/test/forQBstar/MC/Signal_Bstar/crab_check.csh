#!/bin/tcsh

#foreach coup (f-1p0 f-0p5 f-0p1)
foreach coup (f-0p5 f-0p1)
foreach case (500 1000 1500 2000 2500 3000 3500 4000 4500 5000 )

crab status -d Bstar/crab_job_BstarToGJ_M-${case}_${coup}/ >> crab_check.txt

end
end
