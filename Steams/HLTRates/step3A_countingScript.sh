#!/bin/bash

cd CMSSW_14_0_5_patch2/src
eval `scramv1 runtime -sh`
home=$PWD
cd SteamRatesEdmWorkflow/Rates/


#lumiJsons="/afs/cern.ch/work/s/savarghe/public/Run3Rates/json_370293.txt"
lumiJsons="/afs/cern.ch/work/s/savarghe/public/Run3Rates/json_2023D.txt"
#lumiJsons=" /afs/cern.ch/work/s/savarghe/public/Run3Rates/eraGZBlumi.txt"


#Edit config_makeCondorJobsData.py - Input Files Dir
file="config_makeCondorJobsData.py"
line_number=26
new_line="inputFilesDir = \"${home}/SteamRatesEdmWorkflow/Prod/output\""
sed -i "${line_number}s|.*|$new_line|" "$file"

#Edit config_makeCondorJobsData.py - 	CMSSW  Dir
file="config_makeCondorJobsData.py"
line_number=37
new_line="cmsswDir = \"${home}\""
sed -i "${line_number}s|.*|$new_line|" "$file"

#Edit config_makeCondorJobsData.py - LumiJson
file="config_makeCondorJobsData.py"
line_number=42
new_line="json_file = \"${lumiJsons}\""
sed -i "${line_number}s|.*|$new_line|" "$file"


# Once you have made the necessary changes in the config_makeCondorJobsData.py
# (json file, input file directory, cmssw location), create the jobs
python3 config_makeCondorJobsData.py
