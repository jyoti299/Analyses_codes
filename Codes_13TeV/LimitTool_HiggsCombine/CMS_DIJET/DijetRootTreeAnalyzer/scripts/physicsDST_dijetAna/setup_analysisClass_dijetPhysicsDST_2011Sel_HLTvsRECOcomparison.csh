ln -sf $DIJETMACRO/src2011/analysisClass_dijetPhysicsDST_2011Sel_HLTvsRECOcomparison.C $DIJETANA/src/analysisClass.C 
./scripts/make_rootNtupleClass.sh -f /data/santanas/Ntuples/2011/BigNtuplesForTests/RootNtuple-V00-02-08-DATA2011B-AODPlusPhysicsDST-HT-180250_04012012/RootTupleMakerV2_output_DATA_AOD-PhysicsDSTStream2011_2011B-HT-180250_20kevents_0.root -t rootTupleTree/tree
make clean
make
#./main config/RootNtuple-V00-02-08-DATA2011B-AODPlusPhysicsDST-HT-180250_04012012/RootTupleMakerV2_output_DATA_AOD-PhysicsDSTStream2011_2011B-HT-180250_20kevents.txt $DIJETMACRO/config2011/cutTable_dijetPhysicsDST_2011Sel_HLTvsRECOcomparison.txt rootTupleTree/tree analysisClass_dijetPhysicsDST_2011Sel_HLTvsRECOcomparison analysisClass_dijetPhysicsDST_2011Sel_HLTvsRECOcomparison
