import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # TTbar PU200 small BS
        "/store/relval/CMSSW_10_6_0_patch2/RelValTTbar14TeVBSz2p3/GEN-SIM-RECO/PU25ns_106X_upgrade2023_realistic_v3_2023D41PU200-v1/20003/A8405CDA-0657-3C44-82BD-E75A78DDC74F.root"
        )
)

process.analysis = cms.EDAnalyzer(
    'MuonIsolationAnalyzer',
    VertexTag3D  = cms.InputTag("offlinePrimaryVertices"),
    VertexTag4D  = cms.InputTag("offlinePrimaryVertices4D"),
    PileUpTag    = cms.InputTag("addPileupInfo"), 
    muonsTag   = cms.untracked.InputTag("muons", "", "RECO"),
    PFCandidateTag = cms.InputTag("particleFlow", "", "RECO"),
    TrackFastSimTimeValueMapTag = cms.InputTag("trackTimeValueMapProducer", "generalTracksConfigurableFlatResolutionModel", "RECO"),
    TrackFastSimTimeErrValueMapTag = cms.InputTag("trackTimeValueMapProducer", "generalTracksConfigurableFlatResolutionModelResolution", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    mtd5sample = cms.untracked.bool(False),
    #timeResolutions = cms.untracked.vdouble(0.030, 0.040, 0.050, 0.060, 0.070, 0.200),
    timeResolutions = cms.untracked.vdouble(0.030),
    isoConeDR = cms.untracked.double(0.3),
    saveTracks = cms.untracked.bool(True),
    maxDz = cms.untracked.double(0.1),
    maxDxy = cms.untracked.double(9999.), # was 0.02
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

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonIsolation.root"))

process.p = cms.Path(process.analysis)
