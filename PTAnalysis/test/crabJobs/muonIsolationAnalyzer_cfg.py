import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90003/5ADE290E-1F66-134D-806E-CB620FA0EF3E.root' 
        #'/store/mc/PhaseIIMTDTDRAutumn18DR/TTbar_14TeV_TuneCP5_Pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/80003/FFF35381-A845-4841-9364-923285FFCFA6.root'
    ),
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
    mtd5sample = cms.untracked.bool(True),
    timeResolutions = cms.untracked.vdouble(0.035, 0.040, 0.050, 0.060, 0.070, 0.200),
    #timeResolutions = cms.untracked.vdouble(0.035),
    isoConeDR = cms.untracked.double(0.3),
    saveTracks = cms.untracked.bool(True),
    maxDz = cms.untracked.double(0.1),
    minDr = cms.untracked.double(0.0),
    minTrackPt = cms.untracked.double(0.0),
    useVertexClosestToGenZ = cms.untracked.bool(True),
    useVertexClosestToGenZT = cms.untracked.bool(False),
    btlEfficiency = cms.untracked.double(1.0),
    etlEfficiency = cms.untracked.double(1.0)
    #btlEfficiency = cms.untracked.double(0.90),
    #etlEfficiency = cms.untracked.double(0.90)
)

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonIsolation.root"))

process.p = cms.Path(process.analysis)
