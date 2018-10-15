import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'/store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/noPU_93X_upgrade2023_realistic_v2-v1/00000/0047F2E2-92AD-E711-B627-002590A371CA.root'
        '/store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00002/FEF5CD96-78B3-E711-B58F-0CC47AA53D82.root '
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
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    timeResolutions = cms.untracked.vdouble(0.030, 0.050, 0.070),
    isoConeDR = cms.untracked.vdouble(0.2, 0.3, 0.4, 0.5),
    saveTracks = cms.untracked.bool(True),
    maxDz = cms.untracked.double(0.1),
    minDr = cms.untracked.double(0.01)
)

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonIsolation.root"))

process.p = cms.Path(process.analysis)
