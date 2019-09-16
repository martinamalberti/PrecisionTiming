from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'TTbar_muIso'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'muonIsolationAnalyzer_cfg.py'

config.Data.inputDataset = '/RelValTTbar14TeVBSz2p3/CMSSW_10_6_0_patch2-PU25ns_106X_upgrade2023_realistic_v3_2023D41PU200-v1/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 180
config.Data.outLFNDirBase = '/store/user/%s/MTD/MuIsoBS2p3/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'test_TTbar_muIso'

config.Site.storageSite = 'T2_CH_CERN'
