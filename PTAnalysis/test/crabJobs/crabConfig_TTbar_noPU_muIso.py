from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'TTbar_noPU_muIso'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'muonIsolationAnalyzer_mtdTrackReco_cfg.py'
config.JobType.pyCfgParams = ['runMTDReco=1', 'runMTDTrackReco=1', 'runTofPID=1', 'runPUIDMVA=1']

config.Data.inputDataset = '/RelValTTbar_Tauola_14TeV/CMSSW_10_4_0_mtd5-103X_upgrade2023_realistic_v2_2023D35noPU-v1/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 300
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 30
config.Data.outLFNDirBase = '/store/user/%s/MTD/1040mtd5/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'test_TTbar_noPU_muIso'
config.Data.allowNonValidInputDataset = True

config.Site.storageSite = 'T2_CH_CERN'
