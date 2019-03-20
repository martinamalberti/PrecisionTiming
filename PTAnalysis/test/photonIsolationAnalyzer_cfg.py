import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/PhaseIIMTDTDRAutumn18DR/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/280000/01A75567-3B61-9D4E-89B7-D7309FDB0438.root' 
 )
)

process.analysis = cms.EDAnalyzer(
    'PhotonIsolationAnalyzer',
    VertexTag3D  = cms.InputTag("offlinePrimaryVertices"),
    VertexTag4D  = cms.InputTag("offlinePrimaryVertices4D"),
    PileUpTag    = cms.InputTag("addPileupInfo"), 
    barrelPhotonsTag   = cms.untracked.InputTag("gedPhotons"),
    endcapPhotonsTag   = cms.untracked.InputTag("photonsFromMultiCl"),
    PFCandidateTag = cms.InputTag("particleFlow", "", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    TrackTimeValueMapTag = cms.InputTag("trackTimeValueMapProducer", "generalTracksConfigurableFlatResolutionModel", "RECO"), # time at vertex smeared by 30 ps (FastSim)
    TrackTimeErrValueMapTag = cms.InputTag("trackTimeValueMapProducer", "generalTracksConfigurableFlatResolutionModelResolution", "RECO"), # time at vertex smeared by 30 ps (FastSim)  
    #TrackTimeValueMapTag = cms.InputTag("tofPID", "t0", "PAT"),
    #TrackTimeErrValueMapTag = cms.InputTag("tofPID", "sigmat0", "PAT"),
    timeResolutions = cms.untracked.vdouble(0.030),
    isoConeDR = cms.untracked.double(0.3),
    saveTracks = cms.untracked.bool(True),
    maxDz = cms.untracked.double(0.1),
    maxDxy = cms.untracked.double(9999.),
    minDr = cms.untracked.double(0.01),
    btlMinTrackPt = cms.untracked.double(0.0),
    etlMinTrackPt = cms.untracked.double(0.0),
    useVertexClosestToGenZ = cms.untracked.bool(False),
    useVertexClosestToGenZT = cms.untracked.bool(True),
    btlEfficiency = cms.untracked.double(1.0),
    etlEfficiency = cms.untracked.double(1.0)
)

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("photonIsolation.root"))

process.p = cms.Path(process.analysis)
