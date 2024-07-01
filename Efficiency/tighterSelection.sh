#!/bin/bash

######################################
ouput_name="GJ_Def_STEAM_L1DPG"
denQuery="GoodLumi==1 & 2.9<DiMu_mass<3.3 & DiMu_Prob>0.005 & abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 & muTag_pt>8 & muTag_L1_match==1  & muProbe_L1_match==1 & muTag_charge+muProbe_charge==0 & L1muProbe_pt>=5 & L1muTag_pt>=5 & muProbe_pt>=9 & muTag_pt>=9 & L1muTag_quality>=8 & L1muProbe_quality>=8"
numQuery=" "
tagPath=HLT_Mu8_v
probePath=HLT_Mu0_L1DoubleMu_v


######################################
ouput_name="NoID"
denQuery="GoodLumi==1 & 2.9<DiMu_mass<3.3 & DiMu_Prob>0.005 & abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 & muTag_pt>8 & muTag_charge+muProbe_charge==0"
numQuery="none"
tagPath=HLT_Mu8_v
probePath=HLT_Mu0_L1DoubleMu_v









######################################
ouput_name="TightID_BothDen"
denQuery="GoodLumi==1 & 2.9<DiMu_mass<3.3 & DiMu_Prob>0.005 & abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 & muTag_pt>8 & muTag_charge+muProbe_charge==0 & muTagtight==1 & muProbetight==1"
numQuery="none"
tagPath=HLT_Mu8_v
probePath=HLT_Mu0_L1DoubleMu_v


######################################
ouput_name="TightID_BothDen_ProbeL1Match"
denQuery="GoodLumi==1 & 2.9<DiMu_mass<3.3 & DiMu_Prob>0.005 & abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 & muTag_pt>8 & muTag_charge+muProbe_charge==0 & muTagtight==1 & muProbetight==1"
numQuery="muProbe_L1_match==1"
tagPath=HLT_Mu8_v
probePath=HLT_Mu0_L1DoubleMu_v


######################################
ouput_name="TightID_BothDen_ProbeL1Match_ProbepT2"
denQuery="GoodLumi==1 & 2.9<DiMu_mass<3.3 & DiMu_Prob>0.005 & abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 & muTag_pt>8 & muProbe_pt>2 & muTag_charge+muProbe_charge==0 & muTagtight==1 & muProbetight==1"
numQuery="muProbe_L1_match==1"
tagPath=HLT_Mu8_v
probePath=HLT_Mu0_L1DoubleMu_v








######################################
ouput_name="SoftID_BothDen"
denQuery="GoodLumi==1 & 2.9<DiMu_mass<3.3 & DiMu_Prob>0.005 & abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 & muTag_pt>8 & muTag_charge+muProbe_charge==0 & muTagsoft==1 & muProbesoft==1"
numQuery="none"
tagPath=HLT_Mu8_v
probePath=HLT_Mu0_L1DoubleMu_v

######################################
ouput_name="SoftID_BothDen_ProbeL1Match"
denQuery="GoodLumi==1 & 2.9<DiMu_mass<3.3 & DiMu_Prob>0.005 & abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 & muTag_pt>8 & muTag_charge+muProbe_charge==0 & muTagsoft==1 & muProbesoft==1"
numQuery="muProbe_L1_match==1"
tagPath=HLT_Mu8_v
probePath=HLT_Mu0_L1DoubleMu_v

######################################
ouput_name="SoftID_BothDen_ProbeL1Match_ProbepT2"
denQuery="GoodLumi==1 & 2.9<DiMu_mass<3.3 & DiMu_Prob>0.005 & abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 & muTag_pt>8 & muProbe_pt>2 & muTag_charge+muProbe_charge==0 & muTagsoft==1 & muProbesoft==1"
numQuery="muProbe_L1_match==1"
tagPath=HLT_Mu8_v
probePath=HLT_Mu0_L1DoubleMu_v






######################################
ouput_name="DimasSuggestion"
denQuery="GoodLumi==1 & 2.9<DiMu_mass<3.3 & DiMu_Prob>0.005 & abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 & muTag_pt>8 & muProbe_pt>4 & muTag_charge+muProbe_charge==0 & muTagloose==1 & muProbeloose==1"
numQuery="muProbe_L1_match==1"
tagPath=HLT_Mu8_v
probePath=HLT_Mu0_L1DoubleMu_v


######################################
ouput_name="DimasSuggestion_NoL1matching"
denQuery="GoodLumi==1 & 2.9<DiMu_mass<3.3 & DiMu_Prob>0.005 & abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 & muTag_pt>8 & muProbe_pt>4 & muTag_charge+muProbe_charge==0 & muTagloose==1 & muProbeloose==1"
numQuery="none"
tagPath=HLT_Mu8_v
probePath=HLT_Mu0_L1DoubleMu_v


######################################
ouput_name="DimasSuggestion_GlobalMuon"
denQuery="GoodLumi==1 & 2.9<DiMu_mass<3.3 & DiMu_Prob>0.005 & abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 & muTag_pt>8 & muProbe_pt>4 & muTag_charge+muProbe_charge==0 & muTagloose==1 & muProbeloose==1 & muTagGlobal==1 & muProbeGlobal==1"
numQuery="muProbe_L1_match==1"
tagPath=HLT_Mu8_v
probePath=HLT_Mu0_L1DoubleMu_v



#####################################
ouput_name="DimasSuggestion_GlobalMuon_NoL1Matching"
denQuery="GoodLumi==1 & 2.9<DiMu_mass<3.3 & DiMu_Prob>0.005 & abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 & muTag_pt>8 & muProbe_pt>4 & muTag_charge+muProbe_charge==0 & muTagloose==1 & muProbeloose==1 & muTagGlobal==1 & muProbeGlobal==1"
numQuery=""
tagPath=HLT_Mu8_v
probePath=HLT_Mu0_L1DoubleMu_v




# 2022F:
/eos/user/h/hcrottel/pycms_root_37/bin/python TriggerEfficiency_v1.py \
                                    -i /eos/cms/store/group/phys_bphys/trigger/Run2022/TagTuples/Muon_Run2022F_v2_bothNewG.root\
                                    -o $ouput_name\
                                    -r 2022F\
                                    -t $tagPath\
                                    -p $probePath\
                                    --denQ "${denQuery}"\
                                    --numQ "${numQuery}"


# 2023D:
/eos/user/h/hcrottel/pycms_root_37/bin/python TriggerEfficiency_v1.py \
                                    -i /eos/cms/store/group/phys_bphys/trigger/Run2023/TagTuples/Muon_Run2023D_v2_bothNewG.root\
                                    -o $ouput_name\
                                    -r 2023D\
                                    -t $tagPath\
                                    -p $probePath\
                                    --denQ "${denQuery}"\
                                    --numQ "${numQuery}"



#2024C:
/eos/user/h/hcrottel/pycms_root_37/bin/python TriggerEfficiency_v1.py \
                                    -i /eos/cms/store/group/phys_bphys/trigger/Run2024/TagTuples/Muon_Run2024C_v2_bothNewG.root\
                                    -o $ouput_name\
                                    -r 2024C\
                                    -t $tagPath\
                                    -p $probePath\
                                    --denQ "${denQuery}"\
                                    --numQ "${numQuery}"

# 2024D:
/eos/user/h/hcrottel/pycms_root_37/bin/python TriggerEfficiency_v1.py \
                                    -i /eos/cms/store/group/phys_bphys/trigger/Run2024/TagTuples/Muon_Run2024D_v2_bothNewG.root\
                                    -o $ouput_name\
                                    -r 2024D\
                                    -t $tagPath\
                                    -p $probePath\
                                    --denQ "${denQuery}"\
                                    --numQ "${numQuery}"


# 2024E:
/eos/user/h/hcrottel/pycms_root_37/bin/python TriggerEfficiency_v1.py \
                                    -i /eos/cms/store/group/phys_bphys/trigger/Run2024/TagTuples/Muon_Run2024E_v2_bothNewG.root\
                                    -o $ouput_name\
                                    -r 2024E\
                                    -t $tagPath\
                                    -p $probePath\
                                    --denQ "${denQuery}"\
                                    --numQ "${numQuery}"








