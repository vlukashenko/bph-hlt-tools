#!/bin/bash

python3 DimuonHistos.py --input "/eos/user/h/hcrottel/BPHTriggerTuples/ParkingSingleMuon[2-7]/*/*/*/*root"   --outname ParkSinlgeMu2_7_2024B --makePlots
python3 DimuonHistos.py --input "/eos/user/h/hcrottel/BPHTriggerTuples/ParkingDoubleMuonLowMass[2-7]/*/*/*/*root"   --outname ParkDoubleMu2_7_2024B --makePlots

