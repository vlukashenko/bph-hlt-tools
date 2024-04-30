#!/bin/bash

cd CMSSW_14_0_5_patch2/src
eval `scramv1 runtime -sh`
home=$PWD
cd SteamRatesEdmWorkflow/Rates/

./sub_total.jobb
