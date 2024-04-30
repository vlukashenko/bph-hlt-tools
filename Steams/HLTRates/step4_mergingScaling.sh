#!/bin/bash

cd CMSSW_14_0_5_patch2/src
eval `scramv1 runtime -sh`
home=$PWD
cd SteamRatesEdmWorkflow/Rates/


#Edit config_makeCondorJobsData.py - Input Files Dir
file="config_mergeOutputsData.py"
line_number=16
new_line="lumi_in = 1.51"
sed -i "${line_number}s|.*|$new_line|" "$file"

line_number=20
new_line="lumi_target = 2.1"
sed -i "${line_number}s|.*|$new_line|" "$file"

line_number=24
new_line="hlt_ps = 102002"
sed -i "${line_number}s|.*|$new_line|" "$file"

#  Once you modified the config_mergeOutputsData.py as above, run the merging script
python3 config_mergeOutputsData.py
