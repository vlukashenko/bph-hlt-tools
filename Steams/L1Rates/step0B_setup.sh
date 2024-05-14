#!/bin/bash

#menu="https://raw.githubusercontent.com/cms-l1-dpg/L1MenuRun3/master/development/L1Menu_Collisions2024_v1_2_1/L1Menu_Collisions2024_v1_2_1.xml"
menu="https://raw.githubusercontent.com/cms-l1-dpg/L1MenuRun3/master/development/L1Menu_Collisions2024_v1_2_0/L1Menu_Collisions2024_v1_2_0.xml"
filename=$(basename $menu)



cd CMSSW_14_0_4/src
eval `scramv1 runtime -sh`

cd L1MenuTools/rate-estimation

# 3. Translating a menu XML file into C++ code
wget $menu # alternatively: place your custom menu XML here
#wget https://raw.githubusercontent.com/cms-l1-dpg/L1MenuRun3/master/development/L1Menu_Collisions2023_v1_3_0_for2024_v2/L1Menu_Collisions2023_v1_3_0_for2024_v2.xml # alternatively: place your custom menu XML here


#bash configure.sh  L1Menu_Collisions2023_v1_3_0.xml # alternatively: provide your custom menu XML
bash configure.sh $filename # alternatively: provide your custom menu XML




# 4. Compile the rate estimation framework with your custom menulib.* files
mkdir -p objs/include
make -j 8
