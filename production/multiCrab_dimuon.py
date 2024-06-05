#!/usr/bin/env python
"""
This is a small script that does the equivalent of multicrab.
source /cvmfs/cms.cern.ch/common/crab-setup.sh
"""
import os
from optparse import OptionParser

import CRABClient
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
#from httplib import HTTPException #python2
from http.client import HTTPException #python3
#from CRABClient.UserUtilities import config


_ProductionTag = '_First2024Check_April24_TEST'
_ProductionTag = '_MuonJsonTest'
_ProductionTag = "HLTStudy"
def getOptions():
    """
    Parse and return the arguments provided by the user.
    """
    usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS]"
             "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
             "\nUse multicrab -h for help")

    parser = OptionParser(usage=usage)

    parser.add_option('-c', '--crabCmd',
                      dest = 'crabCmd',
                      default = '',
                      help = "crab command",
                      metavar = 'CMD')

    parser.add_option('-w', '--workArea',
                      dest = 'workArea',
                      default = '',
                      help = "work area directory (only if CMD != 'submit')",
                      metavar = 'WAD')

    parser.add_option('-o', '--crabCmdOpts',
                      dest = 'crabCmdOpts',
                      default = '',
                      help = "options for crab command CMD",
                      metavar = 'OPTS')

    (options, arguments) = parser.parse_args()

    if arguments:
        parser.error("Found positional argument(s): %s." % (arguments))
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if options.crabCmd != 'submit':
        if not options.workArea:
            parser.error("(-w WAR, --workArea=WAR) option not provided.")
        if not os.path.isdir(options.workArea):
            parser.error("'%s' is not a valid directory." % (options.workArea))

    return options


def main():

    options = getOptions()
    
    # The submit command needs special treatment.
    if options.crabCmd == 'submit':
        
        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        from CRABClient.UserUtilities import config
        config = config()
        
        #jobname = 'Bsmumuphi_miniAOD_LowMass1_2023Cv0'
        
        config.General.requestName = None
        #config.General.workArea = 'TriggerEfficiencies_Try2'
        config.General.workArea = 'BPHTriggerTuples'
        config.General.transferOutputs = True
        config.General.transferLogs = False

        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = '../test/mumuRootupler.py'
	    #config.JobType.allowUndistributedCMSSW = True
        
        config.Data.inputDataset = None
        config.Data.inputDBS = 'global'
        #config.Data.splitting = 'FileBased'
        config.Data.splitting = "LumiBased"
        config.Data.unitsPerJob = 5
        #config.Data.totalUnits  = 500
        #config.Data.lumiMask = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/DCSOnly_JSONS/dailyDCSOnlyJSON/Collisions24_13p6TeV_378981_380403_DCSOnly_TkPx.json'
        #config.Data.runRange = '379765-379769'

        config.Data.outLFNDirBase = '/store/user/hcrottel/'+str(config.General.workArea)
        #config.Data.outLFNDirBase = '/store/group/phys_bphys/trigger/Run2022/'#+str(config.General.workArea)
        config.Data.publication = False
        config.Data.outputDatasetTag = None
        config.Site.storageSite = 'T3_CH_CERNBOX'
        #config.Site.whitelist = ['T2_US*']
        #config.Data.ignoreLocality = True
        #config.Site.storageSite = None # Choose your site.
        
        #--------------------------------------------------------
        # Will submit one task for each of these input datasets.
        inputDatasets = [
            # "/Muon0/Run2023C-PromptReco-v1/MINIAOD",
            # "/Muon1/Run2023C-PromptReco-v1/MINIAOD",

            # "/Muon0/Run2023C-PromptReco-v2/MINIAOD",
            # "/Muon1/Run2023C-PromptReco-v2/MINIAOD",

            # "/Muon0/Run2023C-PromptReco-v3/MINIAOD",
            # "/Muon1/Run2023C-PromptReco-v3/MINIAOD",

            # "/Muon0/Run2023C-PromptReco-v4/MINIAOD",
            # "/Muon1/Run2023C-PromptReco-v4/MINIAOD",

            # "/Muon0/Run2023D-PromptReco-v1/MINIAOD",
            # "/Muon1/Run2023D-PromptReco-v1/MINIAOD",

            # "/Muon0/Run2023D-PromptReco-v2/MINIAOD",
            # "/Muon1/Run2023D-PromptReco-v2/MINIAOD",

            #"/Muon/Run2022F-PromptReco-v1/MINIAOD",

            #"/ParkingDoubleElectronLowMass/Run2023C-PromptReco-v1/MINIAOD",
            #"/ParkingDoubleElectronLowMass/Run2023C-PromptReco-v2/MINIAOD",
            #"/ParkingDoubleElectronLowMass/Run2023C-PromptReco-v3/MINIAOD",
            #"/ParkingDoubleElectronLowMass/Run2023C-PromptReco-v4/MINIAOD",

            #'/ParkingDoubleMuonLowMass0/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass1/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass2/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass3/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass4/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass5/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass6/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass7/Run2023C-PromptReco-v4/MINIAOD',

            #'/ParkingDoubleMuonLowMass0/Run2023D-PromptReco-v1/MINIAOD',
            #'/ParkingDoubleMuonLowMass2/Run2023D-PromptReco-v1/MINIAOD',
            #'/ParkingDoubleMuonLowMass3/Run2023D-PromptReco-v1/MINIAOD',

            # '/ParkingDoubleMuonLowMass0/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingDoubleMuonLowMass1/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingDoubleMuonLowMass2/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingDoubleMuonLowMass3/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingDoubleMuonLowMass4/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingDoubleMuonLowMass5/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingDoubleMuonLowMass6/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingDoubleMuonLowMass6/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingDoubleMuonLowMass0/Run2024C-PromptReco-v1/MINIAOD',
            # '/ParkingDoubleMuonLowMass1/Run2024C-PromptReco-v1/MINIAOD',
            # '/ParkingDoubleMuonLowMass2/Run2024C-PromptReco-v1/MINIAOD',
            # '/ParkingDoubleMuonLowMass3/Run2024C-PromptReco-v1/MINIAOD',
        

            # '/ParkingSingleMuon0/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingSingleMuon1/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingSingleMuon2/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingSingleMuon3/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingSingleMuon4/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingSingleMuon5/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingSingleMuon6/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingSingleMuon7/Run2024B-PromptReco-v1/MINIAOD',
            # '/ParkingSingleMuon0/Run2024C-PromptReco-v1/MINIAOD',
            # '/ParkingSingleMuon1/Run2024C-PromptReco-v1/MINIAOD',

            # "/Muon0/Run2024C-PromptReco-v1/MINIAOD",
            # "/Muon1/Run2024C-PromptReco-v1/MINIAOD",

            "/Muon0/Run2024D-PromptReco-v1/MINIAOD",
            "/Muon1/Run2024D-PromptReco-v1/MINIAOD",

            #"/Muon0/Run2024E-PromptReco-v1/MINIAOD",
            #"/Muon1/Run2024E-PromptReco-v1/MINIAOD",


            # "/BTagMu/Run2024C-PromptReco-v1/MINIAOD",
            # "/DisplacedJet/Run2024C-PromptReco-v1/MINIAOD",
            # "/EGamma0/Run2024C-PromptReco-v1/MINIAOD",
            # "/JetMET0/Run2024C-PromptReco-v1/MINIAOD",
            # "/Muon0/Run2024C-PromptReco-v1/MINIAOD",
            # "/MuonEG/Run2024C-PromptReco-v1/MINIAOD",
            # "/ParkingHH/Run2024C-PromptReco-v1/MINIAOD",
            # "/ParkingLLP/Run2024C-PromptReco-v1/MINIAOD",
            # "/ParkingVBF0/Run2024C-PromptReco-v1/MINIAOD",
            # "/Tau/Run2024C-PromptReco-v1/MINIAOD",


            # "/Muon0/Run2023D-PromptReco-v1/MINIAOD",
            # "/Muon1/Run2023D-PromptReco-v1/MINIAOD",
            # "/Muon0/Run2023D-PromptReco-v2/MINIAOD",
            # "/Muon1/Run2023D-PromptReco-v2/MINIAOD",

            #"/Muon/Run2022F-PromptReco-v1/MINIAOD",
            #"/ParkingDoubleMuonLowMass0/Run2023D-PromptReco-v1/MINIAOD",
            #"/ParkingDoubleMuonLowMass0/Run2023D-PromptReco-v2/MINIAOD",
        ]
        
        for inDS in inputDatasets:
            # inDS is of the form /A/B/C. Since B is unique for each inDS, use this in the CRAB request name.
            #config.General.requestName = inDS.split('/')[1]+'-'+inDS.split('/')[2]
            config.General.requestName = inDS.split('/')[1]+'-'+inDS.split('/')[2] + _ProductionTag#
            config.Data.inputDataset = inDS
            config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)
            # Submit.
            try:
                print ( "Submitting for input dataset %s" % (inDS) )
                crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print ( "Submission for input dataset %s failed: %s" % (inDS, hte.headers) )
            except ClientException as cle:
                print ( "Submission for input dataset %s failed: %s" % (inDS, cle) )
                
                # All other commands can be simply executed.
    elif options.workArea:

        for dir in os.listdir(options.workArea):
            projDir = os.path.join(options.workArea, dir)
            if not os.path.isdir(projDir):
                continue
            # Execute the crab command.
            msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
            print ( "-"*len(msg) )
            print ( msg )
            print ( "-"*len(msg) )
            try:
                crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print ( "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers) )
            except ClientException as cle:
                print ( "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle) )


if __name__ == '__main__':
    main()
