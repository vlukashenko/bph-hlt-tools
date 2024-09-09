from messages import info, warning, error
from helpers import getPrescaleIdx_per_run, getL1RateInfo, getL1TotalRateInfo, getRunInfo, getLatestRun, getHLTRateInfo, getLumiSecLumi, isGoodRun, getActiveComponents, isCollisions, isLuminosity, isInfo 
from generic_helpers import flattenList, removeDuplicates
from plotting import reformat_for_plotting
import matplotlib.pyplot as plt
import mplhep as hep
import json
import yaml
import os
hep.style.use("CMS")

def test_prescale(omsapi, runs, output_path):
    info("Running prescales test")
    prescaleIdx_for_each_run = []
    for r in runs:
        if not isGoodRun(omsapi, r): continue
        run_recorded_lumi, run_initial_prescale_index, run_components = getRunInfo(omsapi, r)
        if run_recorded_lumi == None:
            continue
        buffer = getPrescaleIdx_per_run(omsapi, r)   
        if buffer == None:
            continue
        prescaleIdx_for_each_run.append({'run_id': r, 'prescale_idx' : buffer})

    fig, ax1 = plt.subplots(figsize=(150,10))
    plot_inputs = reformat_for_plotting(prescaleIdx_for_each_run)

    ax1.plot(range(len(plot_inputs['lumisection'])), plot_inputs['prescale_indices'], 'bo')
    ax2 = ax1.twinx()
    ax2.plot(range(len(plot_inputs['lumisection'])), plot_inputs['recorded_lumi_per_lumisection'], 'ro')

    ax1.set_xticks(range(len(plot_inputs['lumisection'])))
    ax1.set_xticklabels(plot_inputs['lumisection'])
    ax1.tick_params(axis='x', rotation=90)

    plt.xlabel('lumisection')
    ax1.set_ylabel('prescale_indices')
    ax2.set_ylabel('Collected Lumi invpb')
    ax2.spines['right'].set_position(('outward', 0))
    ax2.tick_params(axis='y', colors='r')
    plt.title('Runs {} - {}'.format(runs[0], runs[-1]))
    plot_name = lambda extension : output_path + '/test_prescales.'+extension
    hep.cms.label(rlabel="")
    plt.savefig(plot_name('png'))
    plt.savefig(plot_name('pdf'))

    info('Saving output as ' + plot_name('png') + ' and .pdf')

           
       
def test_l1rate(omsapi, runs, triggers, output_path):
    info("Running total l1rate versus instanteneous lumi test")
    l1_rates = []
    collected_lumi = {}
    harvested_data_per_run = {}
    for r in runs:
        if not isGoodRun(omsapi,r): continue
        
        harvested_data = getL1RateInfo(omsapi, r, triggers)
        if not harvested_data:
            info("Could not find triggers requested active in Run {}".format(r))
            continue
        lumisections = [ trig["lumisec"] for trig in harvested_data.values()] 
        lumisections = flattenList(lumisections)
        lumisections = removeDuplicates(lumisections)
        for lumisec in lumisections:
              name_lumisec = str(r)+"_"+str(lumisec)
              if r not in collected_lumi.keys():
                  collected_lumi[r] = []
              collected_lumi[r].append(getLumiSecLumi(omsapi, name_lumisec))
        harvested_data_per_run[r] = { "data" : harvested_data, "lumisecinfo" : collected_lumi[r]}
        output_json_name = output_path+"/"+str(r)+"_harvested_data.json"

        with open(output_json_name, "w") as output_file:
            info('Saving data harvested for run {} in {}'.format(str(r), output_json_name))
            json.dump(harvested_data_per_run[r], output_file)
               
    
       
        for trig in harvested_data.keys():
            assert len(harvested_data[trig]["post_rate"]) == len(collected_lumi[r]), error("Length of post_rate list and collected_lumi list for trigger {} is not the same: {} vs {}".format(trig, len(harvested_data[trig]["post_rate"]), len(collected_lumi[r])))
            fig, ax1 = plt.subplots(figsize=(15,10))
 
            ax1.plot( collected_lumi[r], harvested_data[trig]["post_rate"], 'bo')
           
            plt.xlabel('L, [invpb]')
            ax1.set_ylabel('L1 Total Rate post scale')
            plt.text(0.1, max(harvested_data[trig]["post_rate"])-0.1*max(harvested_data[trig]["post_rate"]), '{} \n Run {}'.format(trig, r), fontsize=14)
    
            plot_name = lambda extension : output_path + '/test_'+trig+'_r_'+str(r)+'_rate_versus_lumi_per_lumisec.'+extension
            #hep.style.use({"font.size":12})

            hep.cms.label('Preliminary', data=True, lumi = 0.001*int(sum(collected_lumi[r])), year  = "2024")
            plt.savefig(plot_name('png'))
            plt.savefig(plot_name('pdf'))

            info('Saving output as ' + plot_name('png') + ' and .pdf')
            

       
def test_l1totalrate(omsapi, runs, output_path):
    info("Running total l1rate versus instanteneous lumi test")
    l1_rates = []
    collected_luminosity = []
    for r in runs:
        if not isGoodRun(omsapi,r): continue 
        l1_rate, l1_physics, l1_calibration, l1_random = getL1TotalRateInfo(omsapi, r)
        collected_lumi, _, _ = getRunInfo(omsapi, r)
        info("Total Lumi collected in run {} is {} invpb".format(r, collected_lumi))
        info("Total L1 rate for run {} is {} Hz".format(r, l1_rate))
        info("Physics L1 rate for run {} is {} Hz".format(r, l1_physics)) 
        collected_luminosity.append(collected_lumi)
        l1_rates.append(l1_rate)
    
    fig, ax1 = plt.subplots(figsize=(150,10))
  
    ax1.plot(collected_luminosity, l1_rates, 'bo')

    plt.xlabel('L, [invpb]')
    ax1.set_ylabel('L1 Total Rate')
    plt.title('Runs {} - {}'.format(runs[0], runs[-1]))
    plot_name = lambda extension : output_path + '/test_l1totalrate_vs_collected_lumi.'+extension
    hep.cms.label(rlabel="")
    plt.savefig(plot_name('png'))
    plt.savefig(plot_name('pdf'))

    info('Saving output as ' + plot_name('png') + ' and .pdf')


def test_hltrate(omsapi, triggers, runs, output_path):
    info("Running total l1rate test")
    l1rates = []
    with open(triggers, 'r') as file:
        test_trigger_paths = list(json.load(file).keys())
    for r in runs:
        if not isGoodRun(omsapi, r): continue        
        for t in test_trigger_paths:
            hlt_rate, hlt_physics, hlt_calibration, hlt_random = getHLTRateInfo(omsapi, t, r)
            info("Total HLT rate for trigger {} in run {} is {} Hz".format(t, r, hlt_rate))
            info("Physics HLT rate for trigger {} in run {} is {} Hz".format(t, r, hlt_physics))
   
def test_lumi():
    print("PLACEHOLDER")

def test_bad_runs(omsapi, runs, output_path):
    info("Look for the bad runs. The isGoodRun flag is identified by : run in collisions mode, recorded_lumi > 0, info exists.")
    with open('configuration.yaml', 'r') as f:
        config = yaml.safe_load(f)

    save_info = {} #[isGoodRun, components == config["components"] 
    runs = isCollisions(omsapi, runs) #choose only runs with pp collisions
    runs = isLuminosity(omsapi, runs) #choose runs with good lumi
    runs = isInfo(omsapi, runs) #choose runs with info present
    for r in runs:
        components = getActiveComponents(omsapi,r)
        #daqstate = getDAQState(omsapin, r)
        if isGoodRun(omsapi, r, verbose=True) and components == config['components']: continue
        else:
           info("Run {}".format(r) + " is " + str(isGoodRun(omsapi, r)) + " with " + ",".join(map(str, components)) + " components")
           save_info[r] = [isGoodRun, list(set(config['components']) - set(components))]
        
    with open(output_path+'/bad_runs.yaml', 'w') as f:
        yaml.dump(save_info, f)
        info("Saved bad runs info into: " +  output_path+'/bad_runs.yaml')

    
    info("Found {}/{} bad runs".format(len(save_info.keys()), len(runs)))
     
