import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'/store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/noPU_93X_upgrade2023_realistic_v2-v1/00000/0047F2E2-92AD-E711-B627-002590A371CA.root'
        
        #'/store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00002/FEF5CD96-78B3-E711-B58F-0CC47AA53D82.root',
        #'/store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/0051AD4B-31AF-E711-8324-7845C4FC3683.root',
        #'/store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/005E34EC-46B0-E711-99DE-0025905B8576.root',
        #'/store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/8ACC70B3-4AB0-E711-B53C-003048FFD71C.root',         
        #'/store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/8A430CA9-38B0-E711-A1DC-3417EBE2F0DF.root', 

        '/store/mc/PhaseIITDRFall17DR/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/0C84DC5A-01BC-E711-8786-48D539F38892.root',
        #'/store/mc/PhaseIITDRFall17DR/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/30003/A6B64058-09BD-E711-AE05-90B11CBCFFA9.root',
        #'/store/mc/PhaseIITDRFall17DR/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/30003/A8EA1BB6-FEBC-E711-974A-0025905C53B2.root' 
 
        #'/store/mc/PhaseIISpr18AODMiniAOD/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/AODSIM/PU200_93X_upgrade2023_realistic_v5-v1/90000/FEA5A761-EE47-E811-9C32-0CC47A4C8E38.root' 
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
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    timeResolutions = cms.untracked.vdouble(0.030, 0.040, 0.050, 0.060, 0.070),
    isoConeDR = cms.untracked.vdouble(0.2, 0.3, 0.4, 0.5),
    saveTracks = cms.untracked.bool(True),
    maxDz = cms.untracked.double(0.1),
    minDr = cms.untracked.double(0.0),
    useVertexClosestToGen = cms.untracked.bool(True)
)

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonIsolation.root"))

process.p = cms.Path(process.analysis)
