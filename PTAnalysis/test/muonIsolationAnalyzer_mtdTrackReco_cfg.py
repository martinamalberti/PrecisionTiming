# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --filein dbs:/VHToGG_M125_14TeV_amcatnloFXFX_madspin_pythia8/PhaseIIMTDTDRAutumn18DR-PU200_103X_upgrade2023_realistic_v2-v2/FEVT --fileout file:HIG-PhaseIIMTDTDRAutumn18MiniAOD-00006.root --mc --eventcontent MINIAODSIM --runUnscheduled --datatier MINIAODSIM --conditions 103X_upgrade2023_realistic_v2 --step PAT --nThreads 2 --geometry Extended2023D35 --era Phase2C4_timing_layer_bar --python_filename HIG-PhaseIIMTDTDRAutumn18MiniAOD-00006_1_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring
import subprocess
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('datasets',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "Input dataset(s)")
options.register('output',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "output file name")
options.register('pattern',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "pattern of file names to be processed")
options.register('debug',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Print debug messages")
options.register('useparent',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Load data from parent datasets")
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
options.register('runPUIDMVA',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run MTD Reco")
options.register('runPAT',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "RunPAT")
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
#    numberOfThreads=cms.untracked.uint32(1),
#    numberOfStreams=cms.untracked.uint32(0),
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

files = []
secondary_files = []
for dataset in options.datasets:
    print('>> Creating list of files from: \n'+dataset)
    for instance in ['global', 'phys03']:
        query = "-query='file dataset="+dataset+" instance=prod/"+instance+"'"
        if options.debug:
            print(query)
        lsCmd = subprocess.Popen(['dasgoclient '+query+' -limit=0'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        str_files, err = lsCmd.communicate()
        files.extend(['root://cms-xrd-global.cern.ch/'+ifile for ifile in str_files.split("\n")])
        files = [k for k in files if options.pattern in k]
        files.pop()
        if options.useparent:
            print('>> Creating list of secondary files from: \n'+dataset)
            for file in files:
                query = "-query='parent file="+file[len('root://cms-xrd-global.cern.ch/'):]+" instance=prod/"+instance+"'"
                if options.debug:
                    print(query)
                lsCmd = subprocess.Popen(['dasgoclient '+query+' -limit=0'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                str_files, err = lsCmd.communicate()
                secondary_files.extend(['root://cms-xrd-global.cern.ch/'+ifile for ifile in str_files.split("\n")])
                secondary_files.pop()
        
if options.debug:
    for ifile in files:
        print(ifile)

# Input source
process.source = cms.Source(
    "PoolSource",
    #fileNames = cms.untracked.vstring(files),
    #secondaryFileNames = cms.untracked.vstring(secondary_files)
    fileNames = cms.untracked.vstring(
        #'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90003/5ADE290E-1F66-134D-806E-CB620FA0EF3E.root'
        '/store/mc/PhaseIIMTDTDRAutumn18DR/TTbar_14TeV_TuneCP5_Pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/80002/D8DBF74C-A132-0643-9FA4-6153697CF9B8.root'
        ) 
    )
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string(options.output),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

process.reconstruction_step = cms.Path()

if (options.runMTDReco or options.run4DVertex or options.runTofPID or options.runMTDTrackReco):
    process.load("Configuration.StandardSequences.Reconstruction_cff")
    process.pfPileUpIso.PFCandidates = cms.InputTag("particleFlowPtrs")
    process.pfNoPileUpIso.bottomCollection = cms.InputTag("particleFlowPtrs")


if options.runMTDReco:
    print ">>> Run MTD LocalReco"
    process.reconstruction_step += process.mtdClusters
    process.reconstruction_step += process.mtdTrackingRecHits
    process.MINIAODSIMoutput.outputCommands.extend(['keep *_mtdClusters_*_PAT','keep *_mtdTrackingRecHits_*_PAT'])
if options.runMTDTrackReco: 
    print ">>> Run MTD Track Extender with vtx constraint"
    process.trackExtenderWithMTD.UseVertex = cms.bool(True) #run trackExtender using vertex constrain
    if (options.useSimVertex):
        print ">>> Use SimVertex"
        process.trackExtenderWithMTD.UseSimVertex = cms.bool(True) #useSimVertex as the vector constrain
    process.trackExtenderWithMTD.DZCut = 0.3
    process.reconstruction_step += process.trackExtenderWithMTD
    process.MINIAODSIMoutput.outputCommands.extend(['keep *_trackExtenderWithMTD_*_PAT'])

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
    process.tofPID.fixedT0Error = cms.double(0.035) #put a constant 0.035 [ns] error for each track (cannot be used with the current training)
    process.reconstruction_step += process.tofPID
    if options.runPAT:
        process.packedPFCandidates.TimeMap = cms.InputTag("tofPID:t0")
        process.packedPFCandidates.TimeErrorMap = cms.InputTag("tofPID:sigmat0")
    process.MINIAODSIMoutput.outputCommands.extend(['keep *_tofPID_*_PAT'])

if options.runPUIDMVA:
    print ">>> Run PUIDMVA"
    from CommonTools.RecoAlgos.trackPUIDMVAProducer_cfi import trackPUIDMVAProducer
    process.trackPUIDMVA = trackPUIDMVAProducer.clone()
    process.trackPUIDMVA.tracksMTDSrc = cms.InputTag('trackExtenderWithMTD','','PAT')
    process.trackPUIDMVA.btlMatchChi2Src = cms.InputTag('trackExtenderWithMTD', 'btlMatchChi2','PAT')
    process.trackPUIDMVA.btlMatchTimeChi2Src = cms.InputTag('trackExtenderWithMTD', 'btlMatchTimeChi2','PAT')
    process.trackPUIDMVA.etlMatchChi2Src = cms.InputTag('trackExtenderWithMTD', 'etlMatchChi2','PAT')
    process.trackPUIDMVA.etlMatchTimeChi2Src = cms.InputTag('trackExtenderWithMTD', 'etlMatchTimeChi2','PAT')
    process.trackPUIDMVA.mtdTimeSrc = cms.InputTag('trackExtenderWithMTD', 'tmtd','PAT')
    process.trackPUIDMVA.pathLengthSrc = cms.InputTag('trackExtenderWithMTD', 'pathLength','PAT')
    process.trackPUIDMVA.t0TOFPIDSrc = cms.InputTag('tofPID', 't0','PAT')
    process.trackPUIDMVA.sigmat0TOFPIDSrc = cms.InputTag('tofPID', 'sigmat0','PAT')

    process.reconstruction_step += process.trackPUIDMVA
    process.MINIAODSIMoutput.outputCommands.extend(['keep *_trackPUIDMVA_*_PAT'])
# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')

# Path and EndPath definitions
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)

process.Flag_trkPOGFilters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path()
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path()
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
process.Flag_HBHENoiseFilter = cms.Path()
process.Flag_trkPOG_toomanystripclus53X = cms.Path()
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path()
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)


# Load analyzer
process.load('PrecisionTiming.PTAnalysis.muonIsolationAnalyzer_cfi')
muIsoAnalyzer = process.MuonIsolationAnalyzer
process.myanalysis = cms.EndPath(muIsoAnalyzer)

# Output TFile
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('muonIsolation.root')
    )


# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.endjob_step,#process.MINIAODSIMoutput_step,
process.myanalysis)

if options.runPAT:
    process.schedule.associate(process.patTask)
    from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
    associatePatAlgosToolsTask(process)

# customisation of the process.
    
# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 
    
#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 
#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

#Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(2)
#process.options.numberOfStreams=cms.untracked.uint32(0)
