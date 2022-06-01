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
        #'/store/relval/CMSSW_10_4_0_mtd3/RelValTTbar_Tauola_14TeV/GEN-SIM-RECO/103X_upgrade2023_realistic_v2_2023D35noPU_1-v1/20000/8B39B217-65E1-234E-BEA1-9E67280731AE.root'
        #'file:/eos/cms/store/group/phys_egamma/sobhatta/egamma_timing_studies/samples/DYToLL_M-50_TuneCP5_14TeV-pythia8_Phase2HLTTDRWinter20DIGI-PU200_pilot2_110X_mcRun4_realistic_v3-v2_GEN-SIM-DIGI-RAW_2022-03-18_19-23-30/output_10.root'
       '/store/relval/CMSSW_12_4_0_pre3/RelValTTbar_14TeV/MINIAODSIM/PU_123X_mcRun4_realistic_v11_2026D88PU200-v1/2580000/c354eb33-0710-4697-959d-6ae6ffa27946.root'
    ),
    #inputCommands=cms.untracked.vstring(
    #    'keep *',
    #    'drop FTLClusteredmNewDetSetVector_mtdClusters_FTLBarrel_RECO',
    #    'drop FTLClusteredmNewDetSetVector_mtdClusters_FTLEndcap_RECO',
    #    'drop MTDTrackingRecHitedmNewDetSetVector_mtdTrackingRecHits__RECO'
    #    )
)

process.analysis = cms.EDAnalyzer(
    'MuonIsolationAnalyzer',
    VertexTag3D  = cms.InputTag("offlineSlimmedPrimaryVertices"),
    VertexTag4D  = cms.InputTag("offlineSlimmedPrimaryVertices4D"),
    PileUpTag    = cms.InputTag("slimmedAddPileupInfo"), 
    muonsTag   = cms.InputTag("slimmedMuons", "", "RECO"),
    PFCandidateTag = cms.InputTag("packedPFCandidates", "", "RECO"),
    genPartTag = cms.InputTag("packedGenParticles", "", "RECO"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),
    genJetsTag = cms.untracked.InputTag("slimmedGenJets", "", "RECO"),
    timeResolutions = cms.untracked.vdouble(-1.),
    isoConeDR = cms.untracked.double(0.3),
    saveTracks = cms.untracked.bool(True),
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
