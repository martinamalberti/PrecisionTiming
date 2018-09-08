import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'/store/mc/PhaseIITDRFall17DR/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/30000/FEC7F72A-2BBC-E711-972D-A0369F83630C.root'
        #'/store/mc/PhaseIITDRFall17DR/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/30000/82C07609-42BC-E711-86CB-001E6779267C.root'
#'/store/mc/PhaseIITDRFall17DR/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/30000/0403AAD6-CCBC-E711-BC7B-0025904C637C.root'
 '/store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/02527FB1-5CAD-E711-979A-24BE05C63681.root',
'/store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/04EB2BA6-66AD-E711-9003-24BE05C4D8C1.root',
'/store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/0A13405A-5CAD-E711-A5F0-E0071B7A48A0.root' 
 )
)

process.analysis = cms.EDAnalyzer(
    'PhotonIsolationAnalyzer',
    VertexTag3D  = cms.InputTag("offlinePrimaryVertices"),
    VertexTag4D  = cms.InputTag("offlinePrimaryVertices4D"),
    PileUpTag    = cms.InputTag("addPileupInfo"), 
    photonsTag   = cms.untracked.InputTag("gedPhotons"),
    TracksTag    = cms.InputTag("generalTracks"),
    TrackTimeValueMapTag = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
    PFCandidateTag = cms.InputTag("particleFlow", "", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    timeResolutions = cms.untracked.vdouble(0.030, 0.050, 0.070),
    isoConeDR = cms.untracked.vdouble(0.2, 0.3, 0.4, 0.5),
    saveTracks = cms.untracked.bool(True),
    maxDz = cms.untracked.double(0.1),
    minDr = cms.untracked.double(0.01)
)

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("photonIsolation.root"))

process.p = cms.Path(process.analysis)
