import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_9_3_2/RelValZEE_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/10000/F4C7038A-F9AA-E711-8695-4C79BA18144F.root' 
        #'/store/relval/CMSSW_9_3_2/RelValZEE_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/10000/00AD560E-FAAA-E711-A9CC-E0071B7A9800.root' 
        )
)

process.analysis = cms.EDAnalyzer(
    'ElectronIsolationAnalyzer',
    VertexTag3D  = cms.InputTag("offlinePrimaryVertices"),
    VertexTag4D  = cms.InputTag("offlinePrimaryVertices4D"),
    PileUpTag    = cms.InputTag("addPileupInfo"), 
    barrelElectronsTag   = cms.untracked.InputTag("gedGsfElectrons"),
    endcapElectronsTag   = cms.untracked.InputTag("ecalDrivenGsfElectronsFromMultiCl"),
    PFCandidateTag = cms.InputTag("particleFlow", "", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    timeResolutions = cms.untracked.vdouble(0.030, 0.050, 0.070),
    isoConeDR = cms.untracked.vdouble(0.2, 0.3, 0.4, 0.5),
    saveTracks = cms.untracked.bool(False),
    maxDz = cms.untracked.double(0.1),
    minDr = cms.untracked.double(0.01)
)

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("electronIsolation.root"))

process.p = cms.Path(process.analysis)
