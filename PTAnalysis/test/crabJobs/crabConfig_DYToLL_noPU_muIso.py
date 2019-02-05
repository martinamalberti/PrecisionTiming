from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DYToLL_noPU_muIso'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'muonIsolationAnalyzer_cfg.py'

config.Data.inputDataset = '/DYToLL_M-50_14TeV_pythia8/PhaseIIMTDTDRAutumn18DR-NoPU_pilot_103X_upgrade2023_realistic_v2_ext1-v1/FEVT'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 60
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 30
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'test_DYToLL_noPU_muIso'

config.Site.storageSite = 'T2_CH_CERN'

