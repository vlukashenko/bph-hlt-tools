#!/bin/bash

cd CMSSW_14_0_5_patch2/src

wget https://gist.githubusercontent.com/sanuvarghese/d1d9e41116b4bcf355cd0ffdc87ac675/raw/aa0c7820c1eac69583390690a6ff228f276a6435/customizeHLTForL1Skim.py

eval `scramv1 runtime -sh`

menu="/users/hcrottel/2024/SingleMuLowpT/GRunV107/V2"
tag="SingleMu_L1Skims"

hltGetConfiguration $menu --data \
                  --process MYHLT \
                  --type GRun \
                  --prescale 2p0E34+ZeroBias+HLTPhysics \
                  --globaltag 140X_dataRun3_HLT_for2024TSGStudies_v1 \
                  --timing \
                  > hlt${tag}.py
									
#--customise customizeHLTForL1Skim.customizePrescaleSeeds \
