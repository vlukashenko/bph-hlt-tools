#!/bin/bash

cd CMSSW_14_0_5_patch2/src
eval `scramv1 runtime -sh`

#tag="_SingleMuon1_Minus_HLT_Mu0_Barrel_L1HP11_v1"
#tag="_SingleMuon1_Minus_HLT_Mu0_Barrel_L1HP10_v1"
#tag="_SingleMuon1_Minus_HLT_Mu10_Barrel_L1HP11_IP6_v1"
#tag="_SingleMuon1_Minus_HLT_Mu9_Barrel_L1HP10_IP6_v1"
#tag="_SingleMuon1_MinusALL"
#tag="_SingleMuon1"
#tag="_SingleMuon1_SSPath"
#tag="_SameSign"
#tag="Integration_ss_jpsi_low"
#tag="Reference_jpsiLowpt"
#tag="AddL1Seed_jpsiLowpt"
tag="SingleMu_L1Skims"
python3 timing/submit.py hlt${tag}.py --l1menu L1Menu_Collisions2024_v1_2_0_xml --tag _${tag}
