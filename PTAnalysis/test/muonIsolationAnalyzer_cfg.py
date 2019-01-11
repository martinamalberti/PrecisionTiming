import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
    ignoreTotal = cms.untracked.int32(1)
)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_10_4_0_mtd3/RelValDYToLL_M_50_14TeV/GEN-SIM-RECO/103X_upgrade2023_realistic_v2_2023D35noPU_1-v1/20000/FB69EE58-C9C6-7641-A94F-F1ADF8D709C3.root'
        '/store/relval/CMSSW_10_4_0_mtd3/RelValTTbar_Tauola_14TeV/GEN-SIM-RECO/103X_upgrade2023_realistic_v2_2023D35noPU_1-v1/20000/8B39B217-65E1-234E-BEA1-9E67280731AE.root'
    ),
    inputCommands=cms.untracked.vstring(
        'keep *',
        'drop FTLClusteredmNewDetSetVector_mtdClusters_FTLBarrel_RECO',
        'drop FTLClusteredmNewDetSetVector_mtdClusters_FTLEndcap_RECO',
        'drop MTDTrackingRecHitedmNewDetSetVector_mtdTrackingRecHits__RECO'
        )
)

process.analysis = cms.EDAnalyzer(
    'MuonIsolationAnalyzer',
    VertexTag3D  = cms.InputTag("offlinePrimaryVertices"),
    VertexTag4D  = cms.InputTag("offlinePrimaryVertices4D"),
    PileUpTag    = cms.InputTag("addPileupInfo"), 
    muonsTag   = cms.untracked.InputTag("muons", "", "RECO"),
    PFCandidateTag = cms.InputTag("particleFlow", "", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    timeResolutions = cms.untracked.vdouble(0.030, 0.040, 0.050, 0.060, 0.070),
    isoConeDR = cms.untracked.double(0.3),
    saveTracks = cms.untracked.bool(False),
    maxDz = cms.untracked.double(0.1),
    minDr = cms.untracked.double(0.0),
    minTrackPt = cms.untracked.double(0.0),
    useVertexClosestToGen = cms.untracked.bool(True),
    btlEfficiency = cms.untracked.double(1.0),
    etlEfficiency = cms.untracked.double(1.0)
)

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonIsolation.root"))

process.p = cms.Path(process.analysis)
