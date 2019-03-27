import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('Analysis')

options.register('runMTDReco',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run MTD Reco")
options.register('runMTDTrackReco',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run MTD Reco")
options.register('run4DVertex',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run MTD Reco")
options.register('runTofPID',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run MTD Reco")
options.register('useSimVertex',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run MTD Reco")
options.maxEvents = -1
options.parseArguments()

from Configuration.StandardSequences.Eras import eras

process = cms.Process('PAT',eras.Phase2C4_timing_layer_bar)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    numberOfThreads=cms.untracked.uint32(2),
    numberOfStreams=cms.untracked.uint32(0),
    wantSummary = cms.untracked.bool(True)
    )

# import of standard configurations                                                                                                                                                                       
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# import of standard configurations                                                                                                                                                                       
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


#process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

# Global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# Geometry
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')

process.load("Geometry.MTDNumberingBuilder.mtdNumberingGeometry_cfi")
process.load("Geometry.MTDNumberingBuilder.mtdTopology_cfi")
process.load("Geometry.MTDGeometryBuilder.mtdGeometry_cfi")
process.load("Geometry.MTDGeometryBuilder.mtdParameters_cfi")
process.mtdGeometry.applyAlignment = cms.bool(False)


#input
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90003/5ADE290E-1F66-134D-806E-CB620FA0EF3E.root', 
        #'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90003/FF6DA96E-A319-F741-A908-89D93E5B87F1.root',
        #'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90003/FF4620CA-27D8-F347-85DC-E850255AB6BA.root' 
        #'/store/mc/PhaseIIMTDTDRAutumn18DR/TTbar_14TeV_TuneCP5_Pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/80003/FFF35381-A845-4841-9364-923285FFCFA6.root'
    ),
)


# re-run MTD track reco with vertex constrain
process.reconstruction_step = cms.Path()

if (options.runMTDReco or options.run4DVertex or options.runTofPID or options.runMTDTrackReco):
    process.load("Configuration.StandardSequences.Reconstruction_cff")
    process.pfPileUpIso.PFCandidates = cms.InputTag("particleFlowPtrs")
    process.pfNoPileUpIso.bottomCollection = cms.InputTag("particleFlowPtrs")


if options.runMTDReco:
    print ">>> Run MTD LocalReco"
    process.reconstruction_step += process.mtdClusters
    process.reconstruction_step += process.mtdTrackingRecHits


if options.runMTDTrackReco:
    print ">>> Run MTD Track Extender with vtx constraint"
    process.trackExtenderWithMTD.UseVertex = cms.bool(True) #run trackExtender using vertex constrain                                                                                                     
               
    if (options.useSimVertex):
        print ">>> Use SimVertex"
        process.trackExtenderWithMTD.UseSimVertex = cms.bool(True) #useSimVertex as the vector constrain                                                                                                  
               
    process.trackExtenderWithMTD.DZCut = 0.3
    process.reconstruction_step += process.trackExtenderWithMTD

if options.run4DVertex:
    print ">>> Run 4D Vertex"
    process.vertexrereco4dTask = cms.Task(
        process.unsortedOfflinePrimaryVertices4D,
        process.trackWithVertexRefSelectorBeforeSorting4D ,
        process.trackRefsForJetsBeforeSorting4D,
        process.offlinePrimaryVertices4D,
        process.offlinePrimaryVertices4DWithBS,
        process.unsortedOfflinePrimaryVertices4DnoPID ,
        process.trackWithVertexRefSelectorBeforeSorting4DnoPID ,
        process.trackRefsForJetsBeforeSorting4DnoPID ,
        process.offlinePrimaryVertices4DnoPID ,
        process.offlinePrimaryVertices4DnoPIDWithBS,
        process.tofPIDfor4DwithPID
        )
    process.reconstruction_step += process.vertexrereco4dTask

if options.runTofPID:
    print ">>> Run TofPID"
    if (not options.run4DVertex):
        process.tofPID.vtxsSrc = cms.InputTag('offlinePrimaryVertices4DnoPID')
    process.tofPID.fixedT0Error = cms.double(0.035) #put a constant 0.035 [ns] error for each track
    process.reconstruction_step += process.tofPID
    process.packedPFCandidates.TimeMap = cms.InputTag("tofPID:t0")
    process.packedPFCandidates.TimeErrorMap = cms.InputTag("tofPID:sigmat0")

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')

# analysis
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
    mtd5sample = cms.untracked.bool(True),
    #timeResolutions = cms.untracked.vdouble(0.030, 0.040, 0.050, 0.060, 0.070, 0.200),
    timeResolutions = cms.untracked.vdouble(0.035),
    isoConeDR = cms.untracked.double(0.3),
    saveTracks = cms.untracked.bool(True),
    maxDz = cms.untracked.double(0.1),
    minDr = cms.untracked.double(0.0),
    btlMinTrackPt = cms.untracked.double(0.0),
    etlMinTrackPt = cms.untracked.double(0.0),
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

process.schedule = cms.Schedule(process.reconstruction_step)

process.p = cms.Path(process.analysis)
