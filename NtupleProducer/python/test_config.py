import os
if os.path.exists("perfTuple.root"): os.remove("perfTuple.root")

import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar 

import sys
from collections import namedtuple
Jets = namedtuple("Jets", "label tag task fatJet")

inputFile = str(sys.argv[-2])
nEvents = int(sys.argv[-1])
print(f"\nRunning over file: {inputFile}\nNumber of events: {nEvents}\n")


process = cms.Process("RESP", eras.Phase2C17I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:{}'.format(inputFile)),
    inputCommands = cms.untracked.vstring("keep *", 
            "drop l1tPFClusters_*_*_*",
            "drop l1tPFTracks_*_*_*",
            "drop l1tPFCandidates_*_*_*",
            "drop l1tTkPrimaryVertexs_*_*_*")
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

process.centralGen = cms.EDFilter("CandPtrSelector", src = cms.InputTag("genParticlesForMETAllVisible"), cut = cms.string("abs(eta) < 2.4"))
process.barrelGen = cms.EDFilter("CandPtrSelector", src = cms.InputTag("genParticlesForMETAllVisible"), cut = cms.string("abs(eta) < 1.5"))
process.genMetCentralTrue = process.genMetTrue.clone(src = cms.InputTag("centralGen"))
process.extraPFStuff.add(
    process.genParticlesForMETAllVisible,
    process.centralGen,
    process.barrelGen,
    process.genMetCentralTrue
)


""" Get and store gen particles """
process.load("PhysicsTools.NanoAOD.genparticles_cff")
from PhysicsTools.NanoAOD.simpleGenParticleFlatTableProducer_cfi import simpleGenParticleFlatTableProducer
process.l1pfgenTable = simpleGenParticleFlatTableProducer.clone(
    src = cms.InputTag("genParticles"),
    name = cms.string("GenParticles"),
    doc = cms.string("gen particles"),
    externalVariables = cms.PSet(),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table
    variables = cms.PSet(
        pt  = Var("pt", float, precision=8),
        phi = Var("phi", float, precision=8),
        eta = Var("eta", float, precision=8),
        mass = Var("mass", float, precision=8),
        pdgId = Var("pdgId", int, doc="PDG code of the gen particle"),
        genPartIdxMother = Var("?numberOfMothers>0?motherRef(0).key():-1", "int16", doc="index of the mother particle"),
        vz = Var("vz", float, precision=8),
        charge = Var("charge", int, doc="charge id"),
        status = Var("status", int, doc="Particle status. 1=stable, 2=decayed, 3=docayed after longlived particle, 21=unstable, 22=from a longlived particle, 23=unstable and from a longlived particle"),
        statusFlags = (Var(
            "statusFlags().isLastCopyBeforeFSR()                  * 16384 +"
            "statusFlags().isLastCopy()                           * 8192  +"
            "statusFlags().isFirstCopy()                          * 4096  +"
            "statusFlags().fromHardProcessBeforeFSR()             * 2048  +"
            "statusFlags().isDirectHardProcessTauDecayProduct()   * 1024  +"
            "statusFlags().isHardProcessTauDecayProduct()         * 512   +"
            "statusFlags().fromHardProcess()                      * 256   +"
            "statusFlags().isHardProcess()                        * 128   +"
            "statusFlags().isDirectHadronDecayProduct()           * 64    +"
            "statusFlags().isDirectPromptTauDecayProduct()        * 32    +"
            "statusFlags().isDirectTauDecayProduct()              * 16    +"
            "statusFlags().isPromptTauDecayProduct()              * 8     +"
            "statusFlags().isTauDecayProduct()                    * 4     +"
            "statusFlags().isDecayedLeptonHadron()                * 2     +"
            "statusFlags().isPrompt()                             * 1      ",
            "uint16", doc=("gen status flags stored bitwise, bits are: "
                "0 : isPrompt, "
                "1 : isDecayedLeptonHadron, "
                "2 : isTauDecayProduct, "
                "3 : isPromptTauDecayProduct, "
                "4 : isDirectTauDecayProduct, "
                "5 : isDirectPromptTauDecayProduct, "
                "6 : isDirectHadronDecayProduct, "
                "7 : isHardProcess, "
                "8 : fromHardProcess, "
                "9 : isHardProcessTauDecayProduct, "
                "10 : isDirectHardProcessTauDecayProduct, "
                "11 : fromHardProcessBeforeFSR, "
                "12 : isFirstCopy, "
                "13 : isLastCopy, "
                "14 : isLastCopyBeforeFSR, ")
            )),
    )
)


process.l1pfJetTable = cms.EDProducer("L1PFJetTableProducer",
    gen = cms.InputTag("ak4GenJetsNoNu"),
    commonSel = cms.string("pt > 0 && abs(eta) < 5.0"),
    drMax = cms.double(0.2),
    minRecoPtOverGenPt = cms.double(0.1),
    jets = cms.PSet(
        # Gen = cms.InputTag("ak4GenJetsNoNu"),
        Gen_sel = cms.string(f"pt > 0"),
    ),
    moreVariables = cms.PSet(
        nDau = cms.string("numberOfDaughters()"),
    ),
)


process.l1pfFatJetTable = cms.EDProducer("L1PFJetTableProducer",
    gen = cms.InputTag("ak8GenJetsNoNu"),
    commonSel = cms.string("pt > 0 && abs(eta) < 5.0"),
    drMax = cms.double(0.4),
    minRecoPtOverGenPt = cms.double(0.1),
    jets = cms.PSet(
        Gen = cms.InputTag("ak8GenJetsNoNu"),
        Gen_sel = cms.string(f"pt > 0"),
    ),
    moreVariables = cms.PSet(
        nDau = cms.string("numberOfDaughters()"),
    ),
)


process.l1pfcandTable = cms.EDProducer("L1PFCandTableProducer",
    commonSel = cms.string("pt > 0.0 && abs(eta) < 10.0"),
    cands = cms.PSet(
    ),
    moreVariables = cms.PSet(
        puppiWeight = cms.string("puppiWeight"),   # commented out as not a property of gen jets so raises error
        pdgId = cms.string("pdgId"),
        charge = cms.string("charge"),
    ),
)


""" Import gen particles, form gen (ak8) jets, and add to task """
from RecoJets.JetProducers.ak8GenJets_cfi import ak8GenJets
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.extraPFStuff.add(process.genParticlesForJetsNoNu)

# Produce AK8 jets from gen particles
ak8GenJetsNoNu = ak8GenJets.clone( src = "genParticlesForJetsNoNu" )
setattr(process, 'ak8GenJetsNoNu', ak8GenJetsNoNu)
# Define the task and add it to the process
ak8GenJetsNoNuTask = cms.Task(ak8GenJetsNoNu)
setattr(process, 'ak8GenJetsNoNuTask', ak8GenJetsNoNuTask)
process.extraPFStuff.add(process.ak8GenJetsNoNuTask)

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
ak4GenJetsNoNu = ak4GenJets.clone(src="genParticlesForJetsNoNu")
setattr(process, "ak4GenJetsNoNu", ak4GenJetsNoNu)
ak4GenJetsNoNuTask = cms.Task(ak4GenJetsNoNu)
setattr(process, "ak4GenJetsNoNuTask", ak4GenJetsNoNuTask)
process.extraPFStuff.add(process.ak4GenJetsNoNuTask)
# setattr(process.l1pfJetTable.jets, "ak4Gen", cms.InputTag("ak4GenJets"))


""" Form AK8 jets on PUPPI candidates """
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
# Produce AK8 jets from PUPPI candidates
ak8PuppiJets = ak8GenJets.clone(src="l1tLayer2Deregionizer:Puppi")
setattr(process, "ak8PuppiJets", ak8PuppiJets)
# Define the task and add it to the process
ak8PuppiJetsTask = cms.Task(ak8PuppiJets)
setattr(process, "ak8PuppiJetsTask", ak8PuppiJetsTask)
process.extraPFStuff.add(process.ak8PuppiJetsTask)
# Add to jet table
setattr(process.l1pfFatJetTable.jets, "ak8Puppi", cms.InputTag("ak8PuppiJets"))

# Produce AK4 jets from PUPPI candidates
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
ak4PuppiJets = ak4GenJets.clone(src="l1tLayer2Deregionizer:Puppi")
setattr(process, "ak4PuppiJets", ak4PuppiJets)
ak4PuppiJetsTask = cms.Task(ak4PuppiJets)
setattr(process, "ak4PuppiJetsTask", ak4PuppiJetsTask)
process.extraPFStuff.add(process.ak4PuppiJetsTask)
setattr(process.l1pfJetTable.jets, "ak4Puppi", cms.InputTag("ak4PuppiJets"))


""" Save PUPPI candidates to the event """
setattr (process.l1pfcandTable.cands, "puppi", cms.InputTag("l1tLayer2Deregionizer:Puppi"))


"""" ADD JETS TO THE JET TABLE """
def addJets(label, tag, task, fatJet=False):
    process.extraPFStuff.add(task)
    setattr(process.l1pfFatJetTable.jets, label, tag) if fatJet else setattr(process.l1pfJetTable.jets, label, tag)

def addJetConstituents(fatJet=False, N=128):
    for i in range(N): # save a max of N daughters (unfortunately 2D arrays are not yet supported in the NanoAOD output module)
        if fatJet:
            for var in "pt", "eta", "phi", "mass", "pdgId":
                setattr(process.l1pfFatJetTable.moreVariables, "dau%d_%s" % (i, var),    # attribute example dau0_pt
                    cms.string( "? numberOfDaughters() > %d ? daughter(%d).%s : -1"  % (i, i, var) )    # value, example 1st iter:  "? numberOfDaughters() > 0 ? daughter(0).pt : -1"
                    )                                                                                   # failing because number of daughters is not greater than 0 for histojets - num of daughters() returning 0
            setattr(process.l1pfFatJetTable.moreVariables, "dau%d_%s" % (i,"vz"), cms.string("? numberOfDaughters() > %d ? daughter(%d).%s : -1"  % (i,i,"vertex.Z")))    # Not relevant for finding seeds
        else:
            for var in "pt", "eta", "phi", "mass", "pdgId":
                setattr(process.l1pfJetTable.moreVariables, "dau%d_%s" % (i, var),    # attribute example dau0_pt
                    cms.string( "? numberOfDaughters() > %d ? daughter(%d).%s : -1"  % (i, i, var) )    # value, example 1st iter:  "? numberOfDaughters() > 0 ? daughter(0).pt : -1"
                    )                                                                                   # failing because number of daughters is not greater than 0 for histojets - num of daughters() returning 0
            setattr(process.l1pfJetTable.moreVariables, "dau%d_%s" % (i,"vz"), cms.string("? numberOfDaughters() > %d ? daughter(%d).%s : -1"  % (i,i,"vertex.Z")))    # Not relevant for finding seeds


sc4Sim = Jets(label="sc4PuppiSim", tag=cms.InputTag("l1tSC4PFL1Puppi"), task=process.L1TPFJetsTask, fatJet=False)    # Should have False
sc4Emu = Jets(label="sc4PuppiEmu", tag=cms.InputTag("l1tSC4PFL1PuppiEmulator"), task=process.L1TPFJetsEmulationTask, fatJet=False)    # wont have False
sc4SimCorr = Jets(label="sc4PuppiSimCorr", tag=cms.InputTag("l1tSC4PFL1PuppiCorrected"), task=process.L1TPFJetsTask, fatJet=False)    # Should have jet mass
sc4EmuCorr = Jets(label="sc4PuppiEmuCorr", tag=cms.InputTag("l1tSC4PFL1PuppiCorrectedEmulator"), task=process.L1TPFJetsEmulationTask, fatJet=False)    # wont have jet mass

sc8Sim = Jets(label="sc8PuppiSim", tag=cms.InputTag("l1tSC8PFL1Puppi"), task=process.L1TPFJetsTask, fatJet=True)    # Should have jet mass
sc8Emu = Jets(label="sc8PuppiEmu", tag=cms.InputTag("l1tSC8PFL1PuppiEmulator"), task=process.L1TPFJetsEmulationTask, fatJet=True)    # wont have jet mass
sc8SimCorr = Jets(label="sc8PuppiSimCorr", tag=cms.InputTag("l1tSC8PFL1PuppiCorrected"), task=process.L1TPFJetsTask, fatJet=True)    # Should have jet mass
sc8EmuCorr = Jets(label="sc8PuppiEmuCorr", tag=cms.InputTag("l1tSC8PFL1PuppiCorrectedEmulator"), task=process.L1TPFJetsEmulationTask, fatJet=True)    # wont have jet mass


addJets(*sc4Sim)    # seeded cone jets
addJets(*sc4Emu)
addJets(*sc4SimCorr)
addJets(*sc4EmuCorr)

addJets(*sc8Sim)    # seeded cone jets
addJets(*sc8Emu)
addJets(*sc8SimCorr)
addJets(*sc8EmuCorr)

addJetConstituents(fatJet=False, N=32)
addJetConstituents(fatJet=True,  N=32)


""" Handle output of data """
# process.p = cms.Path( process.ntuple + process.l1pfjetTable )
process.p = cms.Path( process.l1pfJetTable + process.l1pfFatJetTable + process.l1pfcandTable + process.l1pfgenTable )
process.p.associate(process.extraPFStuff)
process.TFileService = cms.Service("TFileService", fileName = cms.string("perfTuple.root"))

process.outnano = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string("perfNano.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
)
process.end = cms.EndPath(process.outnano)