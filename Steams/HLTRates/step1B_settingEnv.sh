#!/bin/bash

cd CMSSW_14_0_5_patch2/src
home=$PWD

git clone https://github.com/sanuvarghese/SteamRatesEdmWorkflow.git
eval `scramv1 runtime -sh`
cd SteamRatesEdmWorkflow/Prod/

# If no version is given, takes the last one
menu="/users/hcrottel/2024/SingleMuLowpT/GRunV107/V2"


hltGetConfiguration $menu --full --offline --no-output --data \
									 --process MYHLT \
									 --type GRun \
									 --prescale 2p0E34+ZeroBias+HLTPhysics \
									 --globaltag 140X_dataRun3_HLT_for2024TSGStudies_v1 \
									 --max-events -1 \
									 --l1 L1Menu_Collisions2024_v1_2_0_xml \
										> hlt.py


# Now dump the config file,
edmConfigDump hlt.py > hlt_config.py
