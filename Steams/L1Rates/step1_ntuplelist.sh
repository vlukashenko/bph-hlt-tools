#!/bin/bash

#ntuple_location="/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/caruta/condor/2023EphZB_run370293_13_2_0_pre3-v164_1693563055"
#ntuple_locaton="/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/caruta/condor/2023EphZB_run370293_13_2_0_pre3_NewBMTFQuality_1702498642"
ntuple_location="/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/caruta/condor/EphZB_2023D_run370293_13_3_0_NewBMTFQuality_AXO_1706623249"
ntuple_location="/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/caruta/condor/EphZB_2024C_run379660_14_0_4_1713621759"
ntuple_location="/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/caruta/condor/EphZB_2023D_run370293_14_0_4_1712308611"
ntuple_location="/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/caruta/condor/EphZB_2024C_run379617_menuv110_14_0_4_1713634081"
# Just to add information of the data, extracted from:
#   HowToL1TriggerMenu
filelist_name="379617_PU63"
#filelist_name="370293_PU60LumiLevel"
cd CMSSW_14_0_4/src/
eval `scramv1 runtime -sh`
cd L1MenuTools/rate-estimation/ntuple

#python3 makeFileList.py /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/bundocka/condor/reHcalTP_Nu_11_2_105p20p1_1623921599 > Run3_NuGun_MC_ntuples.list
python3 makeFileList.py $ntuple_location > Run3_${filelist_name}.list
