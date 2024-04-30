#!/bin/bash
cmsrel CMSSW_14_0_5
cd CMSSW_14_0_5/src
eval `scramv1 runtime -sh`

#git cms-init
#add option to restrict single-Path jobs of hltIntegrationTests to a subset of triggers
#https://github.com/cms-sw/cmssw/pull/44007
#git cms-merge-topic -u 44007

scram build -j 8
