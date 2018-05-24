#!/usr/bin/env python

import os, sys, re

jdl_template = """universe = vanilla
Notification = never
Executable = jobs/run_limits_condor.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = stats, Data_and_ResonanceShapes/histo_data_mjj_Run2015D_silverJSON_2445invpb_JEC_Summer15_25nsV6_clean.root, Data_and_ResonanceShapes/ResonanceShapes_gg_13TeV_PU30_Spring15.root, Data_and_ResonanceShapes/ResonanceShapes_qg_13TeV_PU30_Spring15.root, Data_and_ResonanceShapes/ResonanceShapes_qq_13TeV_Phys14Spring15Mix.root, $ENV(HOME)/lib/libBAT.so.5, $ENV(HOME)/lib/libBATmodels.so.3, $ENV(HOME)/lib/libBATmtf.so.0
Output = condor_DUMMY_MASS_DUMMY_FS_$(Cluster)_$(Process).stdout
Error = condor_DUMMY_MASS_DUMMY_FS_$(Cluster)_$(Process).stderr
Log = condor_DUMMY_MASS_DUMMY_FS_$(Cluster)_$(Process).log
Arguments = DUMMY_MASS DUMMY_FS DUMMY_NPES DUMMY_JOB
#+LENGTH="SHORT"
Queue DUMMY_NJOBS
"""

def main():

    argc = len(sys.argv)
    res_type = sys.argv[1]
    masses = range(int(sys.argv[2]),int(sys.argv[3])+int(sys.argv[4]),int(sys.argv[4]))
    npes = 0
    njobs = 1
    if argc>5: npes = int(sys.argv[5])
    if argc>6: njobs = int(sys.argv[6])

    for mass in masses:
        # create jdl file
        jdl_filename = 'run_limits_' + str(mass) + '_' + res_type + '.jdl'
        jdl_file = open(jdl_filename,'w')
        jdl_content = jdl_template
        jdl_content = re.sub('DUMMY_MASS',str(mass),jdl_content)
        jdl_content = re.sub('DUMMY_FS',res_type,jdl_content)
        jdl_content = re.sub('DUMMY_NPES',(str(npes) if npes>0 else ''),jdl_content)
        jdl_content = re.sub('DUMMY_JOB',('$(Process)' if njobs>1 else ''),jdl_content)
        jdl_content = re.sub('DUMMY_NJOBS',str(njobs),jdl_content)
        jdl_file.write(jdl_content)
        jdl_file.close()

        # submit Condor jobs
        os.system('condor_submit ' + jdl_filename)


if __name__ == '__main__':
    main()
