#!/bin/bash

# Not needed for HLTPhysicsData
#L1menu="https://raw.githubusercontent.com/cms-l1-dpg/L1MenuRun3/master/development/L1Menu_Collisions2023_v1_3_0_for2024_v2/L1Menu_Collisions2023_v1_3_0_for2024_v2.xml"

#export SCRAM_ARCH=slc7_amd64_gcc12
#export SCRAM_ARCH=slc7_amd64_gcc900
#export SCRAM_ARHC=el8_amd64_gcc12

cmsrel CMSSW_14_0_5_patch2
cd CMSSW_14_0_5_patch2/src
eval `scramv1 runtime -sh`
git cms-init
scram build -j 4
