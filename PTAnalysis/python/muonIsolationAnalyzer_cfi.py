import FWCore.ParameterSet.Config as cms

# analyzer
MuonIsolationAnalyzer = cms.EDAnalyzer(
    'MuonIsolationAnalyzer',
    VertexTag3D  = cms.InputTag("offlinePrimaryVertices"),
    VertexTag4D  = cms.InputTag("offlinePrimaryVertices4D"),
    PileUpTag    = cms.InputTag("addPileupInfo"), 
    muonsTag   = cms.untracked.InputTag("muons", "", "RECO"),
    PFCandidateTag = cms.InputTag("particleFlow", "", "RECO"),
    TrackTimeValueMapTag = cms.InputTag("tofPID", "t0", "PAT"),
    TrackTimeErrValueMapTag = cms.InputTag("tofPID", "sigmat0", "PAT"),
    TrackPUID3DMVAValueMapTag = cms.InputTag("trackPUIDMVA", "puId3DMVA", "PAT"),
    TrackPUID4DMVAValueMapTag = cms.InputTag("trackPUIDMVA", "puId4DMVA", "PAT"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    mtd5sample = cms.untracked.bool(True),
    #timeResolutions = cms.untracked.vdouble(0.030, 0.040, 0.050, 0.060, 0.070, 0.200),
    timeResolutions = cms.untracked.vdouble(0.035),
    isoConeDR = cms.untracked.double(0.3),
    saveTracks = cms.untracked.bool(True),
    maxDz = cms.untracked.double(0.1),
    minDr = cms.untracked.double(0.0),
    btlMinTrackPt = cms.untracked.double(0.7),
    etlMinTrackPt = cms.untracked.double(0.4),
    useVertexClosestToGenZ = cms.untracked.bool(True),
    useVertexClosestToGenZT = cms.untracked.bool(False),
    btlEfficiency = cms.untracked.double(1.0),
    etlEfficiency = cms.untracked.double(1.0)
    #btlEfficiency = cms.untracked.double(0.90),
    #etlEfficiency = cms.untracked.double(0.90)
)
