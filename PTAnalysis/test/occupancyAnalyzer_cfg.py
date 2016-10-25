import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # Lindsey's file with track timing info
        'file:/tmp/malberti/step3.root'
        
        # relval TTbar
        #'/store/relval/CMSSW_8_1_0_pre12/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v8_2023D1PU200-v1/00000/FCE7E405-BB83-E611-BC35-0CC47A4C8E56.root',
        #'/store/relval/CMSSW_8_1_0_pre12/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v8_2023D1PU200-v1/00000/FC06D10B-8283-E611-AC31-0025905B8564.root',
        #'/store/relval/CMSSW_8_1_0_pre12/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v8_2023D1PU200-v1/00000/FA749F50-1684-E611-984C-0CC47A78A3E8.root',
        #'/store/relval/CMSSW_8_1_0_pre12/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v8_2023D1PU200-v1/00000/F489CC00-BC83-E611-AEC4-0CC47A4D7664.root',
        #'/store/relval/CMSSW_8_1_0_pre12/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v8_2023D1PU200-v1/00000/F47808C7-8983-E611-B613-0CC47A4D75EE.root'


        #relval Zee
        #'/store/relval/CMSSW_8_1_0_pre12/RelValZEE_14/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v8_2023D1PU200-v1/00000/B257DBDC-2391-E611-8B3F-0CC47A4C8E56.root',
        #'/store/relval/CMSSW_8_1_0_pre12/RelValZEE_14/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v8_2023D1PU200-v1/00000/B4FE07CA-7491-E611-9F3D-0CC47A4D7630.root',
        #'/store/relval/CMSSW_8_1_0_pre12/RelValZEE_14/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v8_2023D1PU200-v1/00000/C62934D5-2391-E611-9975-0CC47A7C35F4.root',
        
        
    ),
    inputCommands=cms.untracked.vstring(
        'keep *',
        'drop EcalTimeDigisSorted_mix_EBTimeDigi_HLT',
        'drop EcalTimeDigisSorted_mix_EETimeDigi_HLT',
        'drop l1tHGCFETriggerDigisSorted_hgcalTriggerPrimitiveDigiProducer__HLT'
        )
)

process.analysis = cms.EDAnalyzer(
    'OccupancyAnalyzer',
    VertexTag    = cms.InputTag("offlinePrimaryVertices"),
    #Vertex1DTag  = cms.InputTag("offlinePrimaryVertices1D"),
    #Vertex4DTag  = cms.InputTag("offlinePrimaryVertices4D"),
    PileUpTag    = cms.InputTag("addPileupInfo"), 
    TracksTag    = cms.InputTag("generalTracks"),
    #TrackTimeValueMapTag = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
    PFCandidateTag = cms.InputTag("particleFlow")
)

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string("test_occupancy_RelValZEE.root"))
                                   #fileName = cms.string("test_occupancy_RelValTTbar.root"))
                                   fileName = cms.string("test_occupancy_ZMM_Lindsey.root"))

process.p = cms.Path(process.analysis)
