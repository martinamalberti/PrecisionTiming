from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'HHTo2B2G_phoIso'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'photonIsolationAnalyzer_cfg.py'

config.Data.inputDataset = '/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/PhaseIIMTDTDRAutumn18DR-PU200_103X_upgrade2023_realistic_v2-v1/FEVT'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 60
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 30
config.Data.outLFNDirBase = '/store/user/%s/MTD/1040mtd5/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'test_HHTo2B2G_phoIso'

config.Site.storageSite = 'T2_CH_CERN'
