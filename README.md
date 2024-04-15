# bph-hlt-tools
Ntuples producer + Trigger plots for CMS BPAG

```
cmsrel CMSSW_14_0_5
cd CMSSW_14_0_5/src && cmsenv
git clone https://github.com/horace-cl/bph-hlt-tools.git myAnalyzers/bph-hlt-tools
scram b -j8
```
Run a simple test to produce ntuples
```
cmsRun myAnalyzers/bph-hlt-tools/test/mumuRootupler.py
```
