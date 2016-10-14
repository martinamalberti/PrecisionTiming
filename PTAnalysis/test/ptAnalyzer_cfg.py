import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/tmp/malberti/step3.root'
    ),
    inputCommands=cms.untracked.vstring(
        'keep *',
        'drop EcalTimeDigisSorted_mix_EBTimeDigi_HLT',
        'drop EcalTimeDigisSorted_mix_EETimeDigi_HLT',
        'drop l1tHGCFETriggerDigisSorted_hgcalTriggerPrimitiveDigiProducer__HLT'
        )
)

process.analysis = cms.EDAnalyzer(
    'PTAnalyzer',
    VertexTag    = cms.InputTag("offlinePrimaryVertices"),
    Vertex1DTag  = cms.InputTag("offlinePrimaryVertices1D"),
    Vertex4DTag  = cms.InputTag("offlinePrimaryVertices4D"),
    PileUpTag    = cms.InputTag("addPileupInfo"), 
    #SimHitsTag = cms.InputTag("g4SimHits","FastTimerHits","SIM")
    TracksTag    = cms.InputTag("generalTracks"),
    TrackTimeValueMapTag = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel")
)

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test_pt.root"))

process.p = cms.Path(process.analysis)
