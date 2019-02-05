from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'TTbar_muIso'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'muonIsolationAnalyzer_cfg.py'

config.Data.inputDataset = '/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 60
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 30
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'test_TTbar_muIso'
#config.Data.outputDatasetTag = 'test_TTbar_muIso_vtxClosestZT'

config.Site.storageSite = 'T2_CH_CERN'