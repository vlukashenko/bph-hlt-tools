#!/bin/bash

# 1. Setting up the environment 
cmsrel CMSSW_14_0_4
cd CMSSW_14_0_4/src

# 2. Setting up the MenuTools environment (NOTE: no cmsenv before configuring MenuTools, it creates problems with the virtual environment)
git clone --depth 1 https://github.com/cms-l1-dpg/L1MenuTools.git
#cd L1MenuTools/rate-estimation
