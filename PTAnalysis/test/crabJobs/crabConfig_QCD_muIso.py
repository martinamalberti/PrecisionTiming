from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'QCD_muIso'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'muonIsolationAnalyzer_mtdTrackReco_cfg.py'
config.JobType.pyCfgParams = ['runMTDReco=1', 'runMTDTrackReco=1', 'runTofPID=1', 'runPUIDMVA=1']

config.Data.inputDataset = '/QCD_Pt-15To7000_TuneCP5_Flat_14TeV-pythia8/PhaseIIMTDTDRAutumn18DR-PU200_103X_upgrade2023_realistic_v2-v1/FEVT'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 180
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 20
config.Data.outLFNDirBase = '/store/user/%s/MTD/1040mtd5/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'test_QCD_muIso'
config.Data.allowNonValidInputDataset = True

config.Site.storageSite = 'T2_CH_CERN'

