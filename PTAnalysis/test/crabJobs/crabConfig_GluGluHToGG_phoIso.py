from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'GluGluHToGG_phoIso'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'photonIsolationAnalyzer_cfg.py'

config.Data.inputDataset = '/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v2/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 60
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 30
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'test_GluGluHToGG_phoIso'

config.Site.storageSite = 'T2_CH_CERN'