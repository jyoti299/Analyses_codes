#!/bin/tcsh

setenv pwd $PWD

#foreach InDir (Pt170Mass471dEta10TID  Pt170Mass500dEta18TID  Pt170NoMassCutdEta10TID  Pt170NoMassCutdEta18TID  Pt200Mass500dEta10TID  Pt200NoMassCutdEta10TID Pt170Mass500dEta16TID  Pt170Mass530dEta19TID  Pt170NoMassCutdEta16TID  Pt170NoMassCutdEta19TID  Pt200Mass530dEta16TID  Pt200NoMassCutdEta16TID)
foreach InDir (Pt170Mass561dEta20TID )
#foreach InDir (Pt170Mass560dEta22TID)

#foreach InScript( SignalEfficiencies  fHalf_SignalEfficiencies 	CompilePlotsForLimit)
foreach InScript(fHalf_SignalEfficiencies)

sed -e 's|ABCDEF|'${InDir}'|g'  ${pwd}/${InScript}.C > ${pwd}/${InScript}_tmp.C

cp ${pwd}/${InScript}.C ${pwd}/${InScript}_orig.C
mv ${pwd}/${InScript}_tmp.C ${pwd}/${InScript}.C

root -l -q ${InScript}.C
mv ${pwd}/${InScript}_orig.C ${pwd}/${InScript}.C

end
echo "done for ${InDir}"
end

