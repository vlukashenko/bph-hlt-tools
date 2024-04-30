#!/bin/bash

cd CMSSW_14_0_5/src
eval `scramv1 runtime -sh`

#GT="140X_dataRun3_HLT_for2024TSGStudies_v1"
GT="auto:phase1_2024_realistic"
INPUT_EVENTS=500
#TAG="LowPtLera"
#MENU="/users/valukash/HLT_Mu7_Barrel/V3"
TAG="NewLowPtPathsGRun"
#TAG="ReferenceMenuGRun"
#MENU="/users/hcrottel/2024/SingleMuLowpT/devCMS14HLTv120/V1"
MENU="/dev/CMSSW_14_0_0/GRun/V107"
MENU="/users/hcrottel/2024/SingleMuLowpT/GRunV107/V2"
INPUT_FILE="/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/savarghe/L1Skim/2024_Skim/L1Menu_Collisions2024_v0_0_0_Full/L1_1.root"

hltIntegrationTests \
	$MENU \
  -n $INPUT_EVENTS \
  --input $INPUT_FILE \
	--paths HLT*Barrel*L1HP5*,HLT*Barrel*L1HP8* \
	-x "--globaltag ${GT}" \
  -x "--eras Run3 --l1-emulator uGT --l1  L1Menu_Collisions2024_v1_2_0_xml" \
  -d output_hltTests${TAG} | tee output_hltTests${TAG}.txt
