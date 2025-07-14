import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

import os
if os.path.exists("inputs140X.root"): os.remove("inputs140X.root")

import sys
inputFile = str(sys.argv[-2])
nEvents = int(sys.argv[-1])
print(f"\nRunning over file: {inputFile}\nNumber of events: {nEvents}\n")

process = cms.Process("IN", eras.Phase2C17I13M9)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2026D110Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D110_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_mcRun4_realistic_v3', '')

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")
process.load("L1Trigger.TrackTrigger.ProducerSetup_cff")
process.load("L1Trigger.TrackerDTC.ProducerED_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(f"root://xrootd-cms.infn.it/{inputFile}"),
    inputCommands = cms.untracked.vstring(
        'keep *',
        # 'drop l1tPFJets_*_*_*',
        # 'drop l1tPFTaus_*_*_*',
        # 'drop l1tTrackerMuons_*_*_*',
        # 'drop *_hlt*_*_HLT',
        # 'drop triggerTriggerFilterObjectWithRefs_*_*_HLT'
    ),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents))
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
        numberOfThreads = cms.untracked.uint32(4),
        #numberOfStreams = cms.untracked.uint32(4),
)

process.PFInputsTask = cms.Task(
    process.L1TLayer1TaskInputsTask,
    process.L1THGCalTriggerPrimitivesTask,
    process.TTClustersFromPhase2TrackerDigis,
    process.TTStubsFromPhase2TrackerDigis,
    process.TrackerDTCProducer,
    #process.offlineBeamSpot,
    process.l1tTTTracksFromTrackletEmulation,
    process.l1tTTTracksFromExtendedTrackletEmulation,
    process.TTTrackAssociatorFromPixelDigis,
    process.TTTrackAssociatorFromPixelDigisExtended,
    process.SimL1EmulatorTask,
#    process.l1tTkStubsGmt,

)
process.p = cms.Path(
        process.l1tLayer1 +
        process.l1tLayer2Deregionizer +
        process.l1tLayer2EG
)
process.p.associate(process.PFInputsTask)
process.p.associate(process.SimL1EmulatorTask)

process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("inputs140X.root"),
        outputCommands = cms.untracked.vstring("drop *",
            # --- GEN
            "keep *_genParticles_*_*",
            "keep *_ak4GenJetsNoNu_*_*",
            "keep *_genMetTrue_*_*",
            # --- Track TPs
            "keep *_l1tTTTracksFromTrackletEmulation_*_*",
            "keep *_l1tTTTracksFromExtendedTrackletEmulation_*_*",
            "keep *_TTTrackAssociatorFromPixelDigis_*_*",
            "keep *_TTTrackAssociatorFromPixelDigisExtended_*_*",
            # --- Calo TPs
            "keep *_simEcalEBTriggerPrimitiveDigis_*_*",
            "keep *_simHcalTriggerPrimitiveDigis_*_*",
            "keep *_simCaloStage2Layer1Digis_*_*",
            "keep *_simCaloStage2Digis_*_*",
            # --- Muon TPs
            "keep *_simMuonRPCDigis_*_*",
            "keep *_simMuonGEMPadDigis_*_*",
            "keep *_simMuonGEMPadDigiClusters_*_*",
            "keep *_simMuonGEMPadDigiProducer_*_*",
            "keep *_simDtTriggerPrimitiveDigis_*_*",
            "keep *_simCscTriggerPrimitiveDigis_*_*",
            "keep *_simTwinMuxDigis_*_*",
            "keep *_simBmtfDigis_*_*",
            "keep *_simKBmtfStubs_*_*",
            "keep *_simKBmtfDigis_*_*",
            "keep *_simEmtfDigis_*_*",
            "keep *_simOmtfDigis_*_*",
            "keep *_simGmtCaloSumDigis_*_*",
            "keep *_simGmtStage2Digis_*_*",
            "keep *_simEmtfShowers_*_*",
            "keep *_simGmtShowerDigis_*_*",
            "keep *_simCscTriggerPrimitiveDigisRun3_*_*",
            "keep *_simMuonME0PadDigis_*_*",
            "keep *_me0TriggerDigis_*_*",
            "keep *_simMuonME0PseudoReDigisCoarse_*_*",
            "keep *_me0RecHitsCoarse_*_*",
            "keep *_me0TriggerPseudoDigis_*_*",
            "keep *_me0RecHits_*_*",
            "keep *_me0Segments_*_*",
            "keep *_me0TriggerConvertedPseudoDigis_*_*",
            "keep *_simCscTriggerPrimitiveDigisPhase2_*_*",
            "keep *_simGtExtFakeStage2Digis_*_*",
            "keep *_simGtStage2Digis_*_*",
            "keep *_CalibratedDigis_*_*",
            "keep *_dtTriggerPhase2PrimitiveDigis_*_*",
            # --- HGCal TPs
            "keep l1tHGCalTriggerCellBXVector_l1tHGCalVFEProducer_*_*",
            #"keep l1tHGCalTriggerCellBXVector_l1tHGCalConcentratorProducer_*_*",
            "keep l1tHGCalMulticlusterBXVector_l1tHGCalBackEndLayer2Producer_*_*",
            "keep l1tHGCalTowerBXVector_l1tHGCalTowerProducer_*_*",
            # --- GCT reconstruction
            "keep *_l1tEGammaClusterEmuProducer_*_*",
            "keep *_l1tTowerCalibration_*_*",
            "keep *_l1tCaloJet_*_*",
            "keep *_l1tCaloJetHTT_*_*",
            # "keep *_l1tPhase2L1CaloEGammaEmulator_*_*",
            # "keep *_l1tPhase2CaloPFClusterEmulator_*_*",
            # --- GTT reconstruction
            "keep *_l1tVertexFinder_*_*",
            "keep *_l1tVertexFinderEmulator_*_*",
            "keep *_l1tTrackJets_*_*",
            "keep *_l1tTrackJetsExtended_*_*",
            "keep *_l1tTrackFastJets_*_*",
            "keep *_l1tTrackerEtMiss_*_*",
            "keep *_l1tTrackerHTMiss_*_*",
            "keep *_l1tTrackJetsEmulation_*_*",
            "keep *_l1tTrackJetsExtendedEmulation_*_*",
            "keep *_l1tTrackerEmuEtMiss_*_*",
            "keep *_l1tTrackerEmuHTMiss_*_*",
            "keep *_l1tTrackerEmuHTMissExtended_*_*",
            # --- GMT reconstruction
            "keep *_l1tStubsGmt_*_*",
            "keep *_l1tKMTFMuonsGmt_*_*",
            "keep *_l1tFwdMuonsGmt_*_*",
            "keep *_l1tSAMuonsGmt_*_*",
        ),
        compressionAlgorithm = cms.untracked.string('LZMA'),
        compressionLevel = cms.untracked.int32(4),
        dropMetaData = cms.untracked.string('ALL'),
        fastCloning = cms.untracked.bool(False),
        overrideInputFileSplitLevels = cms.untracked.bool(True),
        eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
        SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p")),
)
process.e = cms.EndPath(process.out)

process.schedule = cms.Schedule([process.p,process.e])

# SUPPRESS GEM ERROR:
# %MSG-w GEMClusterProcessor:   CSCTriggerPrimitivesProducer:simCscTriggerPrimitiveDigisRun3  11-Jul-2025 13:10:29 CEST Run: 1 Event: 16420
# Encountered unphysical GEM pads when making a single cluster, resetting cluster to empty.
# %MSG-w GEMClusterProcessor:   CSCTriggerPrimitivesProducer:simCscTriggerPrimitiveDigis  11-Jul-2025 13:33:20 CEST Run: 1 Event: 15833
# Encountered unphysical GEM pads when making a coincidence cluster, resetting cluster to empty.
process.MessageLogger.cerr.GEMClusterProcessor = cms.untracked.PSet( limit = cms.untracked.int32(0) )

process.out.outputCommands += [ "drop *_l1tHGCalVFEProducer_*_*", ]

if process.L1TMultiJetsTask:
    process.SimL1EmulatorTask.remove(process.L1TMultiJetsTask)
    del process.L1MultiJetProducer
    del process.l1tMultiJetProducerPuppi
    del process.l1tMultiJetProducerPuppiCorrectedEmulator
    del process.L1TMultiJetsTask
    # open("debug_dump_runInputs140X.py", "w").write(process.dumpPython())