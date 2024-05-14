#!/bin/bash
cd /afs/cern.ch/work/h/hcrottel/private/L1Rates_2024/

start=`date +%s`

prescale="https://raw.githubusercontent.com/cms-l1-dpg/L1MenuRun3/master/development/L1Menu_Collisions2024_v1_2_0/PrescaleTable/L1Menu_Collisions2024_v1_2_0.csv"
FILENAME=$(basename $prescale)

L1Prescales="Prescales${FILENAME}"
lumiatble="EphZB2024C_run379617.csv"
filelist_name="Run3_379617_PU63.list"
nevents=$1
name="DoubleMu0er2p0_$1"




wget -O $L1Prescales $prescale
./modify_prescales.sh $L1Prescales L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5 1


cd CMSSW_14_0_4/src/
eval `scramv1 runtime -sh`
cd L1MenuTools/rate-estimation/

cp "../../../../$L1Prescales" menu/
cp "../../../../$lumiatble" menu/




./testMenu2016 \
	 -u menu/$lumiatble \
	 -m menu/$L1Prescales \
   -l ntuple/$filelist_name \
	 -o $name \
   -b 2554 \
	 --lowerPUbound 58 \
	 --upperPUbound 67 \
   --doPlotRate --doPlotEff \
   --maxEvent $nevents \
   --SelectCol 2p0E34 \
   --doPrintPU --doPlotLS 


end=`date +%s`
runtime=$((end-start))

echo $runtime
