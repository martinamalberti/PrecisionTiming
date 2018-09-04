import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/tmp/malberti/step3.root'
        #'/store/mc/PhaseIITDRFall17DR/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/30000/FEC7F72A-2BBC-E711-972D-A0369F83630C.root'
        '/store/mc/PhaseIITDRFall17DR/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/30000/82C07609-42BC-E711-86CB-001E6779267C.root'
    ),
#    inputCommands=cms.untracked.vstring(
#        'keep *',
#        'drop EcalTimeDigisSorted_mix_EBTimeDigi_HLT',
#        'drop EcalTimeDigisSorted_mix_EETimeDigi_HLT',
#        'drop l1tHGCFETriggerDigisSorted_hgcalTriggerPrimitiveDigiProducer__HLT'
 #       )
)

process.analysis = cms.EDAnalyzer(
    'PhotonIsolationAnalyzer',
    #VertexTag    = cms.InputTag("offlinePrimaryVertices"),
    VertexTag    = cms.InputTag("offlinePrimaryVertices4D"),
    PileUpTag    = cms.InputTag("addPileupInfo"), 
    photonsTag   = cms.untracked.InputTag("gedPhotons"),
    TracksTag    = cms.InputTag("generalTracks"),
    TrackTimeValueMapTag = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
    PFCandidateTag = cms.InputTag("particleFlow", "", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    timeResolutions = cms.untracked.vdouble(0.030, 0.050, 0.070),
    isoConeDR = cms.untracked.vdouble(0.2, 0.3, 0.4, 0.5),
    saveTracks = cms.untracked.bool(True)
)

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("photonIsolation.root"))

process.p = cms.Path(process.analysis)
