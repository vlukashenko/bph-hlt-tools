from tools import getOMSAPI, getAppSecret
from helpers import getPrescaleIdx_per_run, getL1RateInfo, getRunInfo, getLatestRun, isGoodRun, getPreviousRun
from generic_helpers import harvestInfoJSON
from messages import info, warning, error
from tests import test_prescale, test_l1totalrate, test_l1rate, test_hltrate, test_lumi, test_bad_runs
import os, sys
import matplotlib.pyplot as plt
import mplhep as hep 
hep.style.use('CMS')

import yaml
import argparse

#General configuration
#First, start OMSAPI connection with kerberos
omsapi = getOMSAPI(getAppSecret())

#Open general configuration
with open('configuration.yaml', 'r') as file:
    configuration = yaml.safe_load(file)

baseline_components = configuration['components']

#run_min = 378238 #first 2024
#run_max = 380712

def controller(args):
    runs = []
    #print(args.era, args.run, args.runmax)
    flag = args.era != None and args.run == -1 and args.runmax == None
    #print(flag)
    if args.run != -1 and args.runmax == None and args.era == None:
        info("Running for the run {}".format(args.run))
        runs = [args.run]
    elif args.runmax != None and args.run == -1 and args.era == None:
        runmax = args.runmax
        if args.runmax == "latest":
            runmax = int(getLatestRun(omsapi))
        info("Running for the runs {}-{}".format(args.runmin, runmax))
        runs = list(range(args.runmin, runmax, 1))
    elif args.run == -1 and args.runmax == None and args.era == None:
        latest_run = getLatestRun(omsapi)
        info("Running for the latest run {}".format(latest_run)) 
        while not isGoodRun(omsapi, latest_run):
           latest_run = getPreviousRun(omsapi, latest_run)
        runs = [int(latest_run)]
    elif args.era != None and args.run == -1 and args.runmax == None:
        runmin, runmax = configuration['eras'][args.era]
        runs = range(runmin, runmax) 
    else:
        info("No run is chosen")

    lumisections = []
    if args.lumisection != None:
        lumisections = [args.lumisection]
        info("Running for the lumisection {}".format(args.lumisection))
    elif args.lumisectionmax != None:
        lumisections = list(range(lumisectionmin, lumisectionmax,1))
        info("Running for the lumisections {}-{}".format(args.lumisectionmin, args.lumisectionmax))
    else:
        info("No lumisection is chosen")

    output_path = args.output_path
    if runs != []:
        output_path = args.output_path + '/'+str(runs[0])+'_'+str(runs[-1])
        os.mkdir(output_path)
    if lumisections != []:
        output_path = args.output_path + '/'+str(lumisections[0])+'_'+str(lumisections[-1])
    os.makedirs(output_path, exist_ok=True)
    if args.test == 'l1rate':
        l1_triggers = harvestInfoJSON(args.triggers, "L1SeedsLogicalExpression")
        formatted_string = ", ".join(l1_triggers)
        info("Found following {} triggers: {}".format(len(l1_triggers), formatted_string))
        test_l1rate(omsapi, runs, l1_triggers, output_path)
    elif args.test == 'l1totalrate':
        test_l1totalrate(omsapi, runs, output_path)

    elif args.test == 'prescale':
        test_prescale(omsapi, runs, output_path)
    elif args.test == 'hltrate':
        test_hltrate(omsapi, args.triggers, runs, output_path)
    elif args.test == 'lumi':
        test_lumi(omsapi, runs, output_path)
    elif args.test == 'runquality':
        test_bad_runs(omsapi, runs, output_path)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='A master code for the BPH trigger studies', formatter_class=argparse.RawTextHelpFormatter)

    # Add arguments
    parser.add_argument('--triggers',  default = 'triggers.json', required = False,
                        help='JSON file with triggers we are interested in. Required.')
    parser.add_argument('--output-path', default = '/eos/user/v/valukash/TriggerContact/OMSRatesNtuple/OMS_ntuplizer/output', 
                        help='Path to the output of the tests')
    parser.add_argument('--test', default = 'l1rate', choices = ['hltrate', 'prescale', 'lumi', 'l1totalrate', 'l1rate', 'runquality'],
            help = 'Type of test executed: \n - l1rate : l1 trigger rate per lumi section (instant lumi); \n - hltrate: placeholder for hltrate study; \n - prescale : prescales per lumisection; \n - lumi: placeholder for lumi plots; \n - l1totalrate: total rate of L1 per collected luminosity; \n - runquality: check the quality of run based on trigger mode, collected lumi and number of active detector components.')
    parser.add_argument('--runmin', type=int, required = False,
                        help = 'First run to check. Range [0, 1000000]')
    parser.add_argument('--runmax', required = False,
                        help = 'Last run to check. If -1 will choose latest available run. Range [-1, 1000000], If "latest" then looks up the latest run automatically.')
    parser.add_argument('--run', default = -1, required = False, 
                        help = 'Specific run to check. By default -1, which will pick the latest available run')
    parser.add_argument('--era', required = False, help = 'Era to check')
    parser.add_argument('--lumisectionmin', required = False,  
                        help = 'First lumisection in range [0 - x]. By default -10000, which means no lumisection choice')
    parser.add_argument('--lumisectionmax', required = False,
                        help = 'Last lumisection in range [0 - x].')
    parser.add_argument('--lumisection',  required = False,
                        help = 'Choose one lumisection in range [0 - x].')
        
    args = parser.parse_args()

    if args.run != -1 and (args.runmin or args.runmax): 
        print("--run and --runmin/--runmax are incompatable. Use help for more info")
        sys.exit(2)
    if args.lumisection and (args.lumisectionmin or args.lumisectionmax): 
        print("--lumisection and --lumisectionmin/--lumisectionmax are incompatable. Use help for more info")
        sys.exit(2)

    controller(args)

