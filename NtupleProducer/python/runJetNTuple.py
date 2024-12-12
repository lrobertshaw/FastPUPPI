import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar
import os

process = cms.Process("RESP", eras.Phase2C17I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/cms/store/cmst3/group/l1tr/FastPUPPI/14_0_X/fpinputs_131X/v9a/VBFHtt_PU200/inputs131X_1.root'),
    inputCommands = cms.untracked.vstring("keep *", 
            "drop l1tPFClusters_*_*_*",
            "drop l1tPFTracks_*_*_*",
            "drop l1tPFCandidates_*_*_*",
            "drop l1tTkPrimaryVertexs_*_*_*",
            "drop l1tKMTFTracks_*_*_*")
)

process.load('Configuration.Geometry.GeometryExtended2026D110Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D110_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi') # needed for HGCAL_noise_fC
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.Configuration.GenMETParticles_cff')
process.load('RecoMET.METProducers.genMetTrue_cfi')

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from RecoMET.METProducers.pfMet_cfi import pfMet

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_mcRun4_realistic_v3', '')

# NOTE: we need this to avoid saving the stubs
process.l1tTrackSelectionProducer.processSimulatedTracks = False

from L1Trigger.L1CaloTrigger.l1tPhase2L1CaloEGammaEmulator_cfi import l1tPhase2L1CaloEGammaEmulator
process.l1tPhase2L1CaloEGammaEmulator = l1tPhase2L1CaloEGammaEmulator.clone()

process.extraPFStuff = cms.Task(
        process.l1tPhase2L1CaloEGammaEmulator,
        process.l1tSAMuonsGmt,
        process.l1tGTTInputProducer,
        process.l1tTrackSelectionProducer,
        process.l1tVertexFinderEmulator,
        process.L1TLayer1TaskInputsTask,
        process.L1TLayer1Task,
        process.L1TLayer2EGTask)

def addJetNTuple(trktype = "extended", nparam = 5):
    # create new jet tupler
    jetColl = "l1tSC4PFL1PuppiExtendedEmulator"
    jetCollCorr = "l1tSC4PFL1PuppiExtendedEmulator"
    if trktype == "baseline":
        jetColl = "l1tSC4PFL1PuppiEmulator"
        jetCollCorr = "l1tSC4PFL1PuppiCorrectedEmulator"

    process.outnano = cms.EDAnalyzer("JetNTuplizer",
        genJets = cms.InputTag("ak4GenJetsNoNu"),
        genParticles = cms.InputTag("genParticles"),
        scPuppiJets = cms.InputTag(jetColl),
        scPuppiJetsCorr = cms.InputTag(jetCollCorr),
        nnTaus = cms.InputTag("l1tNNTauProducerPuppi","L1PFTausNN"),
        genJetsFlavour = cms.InputTag("genFlavourInfo"),
        vtx = cms.InputTag("l1tVertexFinderEmulator","L1VerticesEmulation"),
        multijetIDs = cms.InputTag("l1tMultiJetProducerPuppiCorrectedEmulator", "L1PFMultiJets"),
        electrons = cms.InputTag("l1tLayer2EG","L1CtTkElectron"),
        muons = cms.InputTag("l1tSAMuonsGmt","promptSAMuons"),
    )
    process.endTuple = cms.EndPath(process.outnano)
    outName = "jetTuple_"+trktype+"_"+str(nparam)+".root"
    process.TFileService = cms.Service("TFileService", fileName = cms.string(outName))

# to check available tags:
process.p = cms.Path()
process.p.associate(process.extraPFStuff)
process.p.associate(process.L1TPFJetsExtendedTask)
process.p.associate(process.L1TBJetsTask)
process.p.associate(process.L1TMultiJetsTask)
process.TFileService = cms.Service("TFileService", fileName = cms.string("jetTuple.root"))

def addNNPuppiTaus():
    process.load("L1Trigger.Phase2L1ParticleFlow.L1NNTauProducer_cff")
    process.l1tNNTauProducerPuppi.maxtaus = cms.int32(500)
    process.extraPFStuff.add(process.l1tNNTauProducerPuppi)

def addSeededConeJets():
    process.extraPFStuff.add(process.L1TPFJetsTask)
    process.extraPFStuff.add(process.L1TPFJetsExtendedTask)

def addMultitagging(trktype = "extended"):
    process.load("L1Trigger.Phase2L1ParticleFlow.L1MultiJetProducer_cff")
    if trktype == "extended":
        process.l1tMultiJetProducerPuppiCorrectedEmulator.jets = cms.InputTag("l1tSC4PFL1PuppiExtendedEmulator")
    else:
        process.l1tMultiJetProducerPuppiCorrectedEmulator.jets = cms.InputTag("l1tSC4PFL1PuppiEmulator")
    process.l1tMultiJetProducerPuppiCorrectedEmulator.maxJets = cms.int32(500)
    process.l1tMultiJetProducerPuppiCorrectedEmulator.MultiJetPath = cms.string(os.environ['CMSSW_BASE']+"/src/hls4ml-jettagger/MultiJetBaseline")
    process.extraPFStuff.add(process.L1TMultiJetsTask)

def addGenJetFlavourTable():
    process.load("PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi")
    process.load("PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi")
    process.selectedHadronsAndPartons.partonMode = cms.string("Pythia8")
    process.genFlavourInfo = process.ak4JetFlavourInfos.clone(jets = "ak4GenJetsNoNu")
    process.p += process.selectedHadronsAndPartons
    process.p += process.genFlavourInfo

def goMT(nthreads=2):
    process.options.numberOfThreads = cms.untracked.uint32(nthreads)
    process.options.numberOfStreams = cms.untracked.uint32(0)

if True:
    process.source.fileNames  = [
        'file:/eos/cms/store/cmst3/group/l1tr/FastPUPPI/14_2_X/fpinputs_140X/v0/DYToLL_M10To50_PU200/inputs140X_7.root'
    ] 
    goMT(4)
    trktype = "extended"
    nparam = 5
    addSeededConeJets()
    addMultitagging(trktype = trktype)
    addNNPuppiTaus()
    addGenJetFlavourTable()
    addJetNTuple(trktype = trktype, nparam = nparam)
    if False:
        open("debug_dump_runJetNTuple.py", "w").write(process.dumpPython())
