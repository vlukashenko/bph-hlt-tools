#!/bin/bash

cd CMSSW_14_0_5_patch2/src
eval `scramv1 runtime -sh`
home=$PWD
cd SteamRatesEdmWorkflow/Prod/


cp /tmp/x509up_u134345 /afs/cern.ch/user/h/hcrottel/private/

# Don't forget to change the import inputFileNames 
file="run_steamflow_cfg.py"
line_number=13
new_line="from list_cff_Skim import inputFileNames"
sed -i "${line_number}s/.*/$new_line/" "$file"

line_number=7
new_line="nEvents=200             # number of events to process"
sed -i "${line_number}s/.*/$new_line/" "$file"

cmsRun run_steamflow_cfg.py

line_number=7
new_line="nEvents=-1             # iNew number of events to process"
sed -i "${line_number}s/.*/$new_line/" "$file"

#exit

mkdir -p $home/SteamRatesEdmWorkflow/Prod/output
./cmsCondorData.py run_steamflow_cfg.py\
									$home\
						       		$home/SteamRatesEdmWorkflow/Prod/output\
						       		-n 1 -q longlunch\
									-p /afs/cern.ch/user/h/hcrottel/private/x509up_u134345
