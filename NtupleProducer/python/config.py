import os
if os.path.exists("perfTuple.root"): os.remove("perfTuple.root")

import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar 

import sys
from collections import namedtuple
Jets = namedtuple("Jets", "label tag task")

inputFile = str(sys.argv[-1])    # root://eoscms.cern.ch//eos/user/l/lroberts/P2_Jets/InputData/CMSSW14/TTbar/inputs131X_9.root
wideJets = True
nEvents = -1
print(f"\nRunning over file: {inputFile}\nWide jets: {wideJets}\nNumber of events: {nEvents}\n")

# Handle wide and regular jets
if wideJets == True:
    genJets = "ak8GenJetsNoNu"
    ptCut = 0
elif wideJets == False:
    genJets = "ak4GenJetsNoNu"
    ptCut = 15
else:
    raise Exception("wideJets must be a boolean")

process = cms.Process("RESP", eras.Phase2C17I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:{}'.format(inputFile)),
#     inputCommands = cms.untracked.vstring("keep *", 
#             "drop l1tPFClusters_*_*_*",
#             "drop l1tPFTracks_*_*_*",
#             "drop l1tPFCandidates_*_*_*",
#             "drop l1tTkPrimaryVertexs_*_*_*")
)

process.load('Configuration.Geometry.GeometryExtended2026D95Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D95_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi') # needed for HGCAL_noise_fC
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.Configuration.GenMETParticles_cff')
process.load('RecoMET.METProducers.genMetTrue_cfi')

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
# from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
from RecoMET.METProducers.pfMet_cfi import pfMet

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '131X_mcRun4_realistic_v9', '')

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

if wideJets == True:
    process.load('RecoJets.Configuration.GenJetParticles_cff')
    process.extraPFStuff.add(process.genParticlesForJetsNoNu)

    from RecoJets.JetProducers.ak8GenJets_cfi import ak8GenJets
    ak8GenJetsNoNu = ak8GenJets.clone( src = "genParticlesForJetsNoNu" )
    setattr(process, 'ak8GenJetsNoNu', ak8GenJetsNoNu)

    ak8GenJetsNoNuTask = cms.Task(ak8GenJetsNoNu)
    setattr(process, 'ak8GenJetsNoNuTask', ak8GenJetsNoNuTask)
    process.extraPFStuff.add(process.ak8GenJetsNoNuTask)

    # AK4 jets
    from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
    ak4GenJetsNoNu = ak4GenJets.clone( src = "genParticlesForJetsNoNu" )
    setattr(process, 'ak4GenJetsNoNu', ak4GenJetsNoNu)
    
    ak4GenJetsNoNuTask = cms.Task(ak4GenJetsNoNu)
    setattr(process, 'ak4GenJetsNoNuTask', ak4GenJetsNoNuTask)
    # process.extraPFStuff.add(process.ak4GenJetsNoNuTask)

    # AK8 ON PUPPI CANDS
    from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
    ak8PuppiJets = ak8GenJets.clone(src="l1tLayer2Deregionizer:Puppi")
    setattr(process, "ak8PuppiJets", ak8PuppiJets)

    ak8PuppiJetsTask = cms.Task(ak8PuppiJets) #l1tLayer2Deregionizer, 
    setattr(process, "ak8PuppiJetsTask", ak8PuppiJetsTask)
    #process.extraPFStuff.add(process.ak8PuppiJetsTask)

    # AK4 ON PUPPI CANDS
    ak4PuppiJets = ak4GenJets.clone(src="l1tLayer2Deregionizer:Puppi")
    setattr(process, "ak4PuppiJets", ak4PuppiJets)

    ak4PuppiJetsTask = cms.Task(ak4PuppiJets) #l1tLayer2Deregionizer, 
    setattr(process, "ak4PuppiJetsTask", ak4PuppiJetsTask)

process.ntuple = cms.EDAnalyzer("ResponseNTuplizer",
    genJets = cms.InputTag(genJets),
    genParticles = cms.InputTag("genParticles"),
    isParticleGun = cms.bool(False),
    writeExtraInfo = cms.bool(False),
    doRandom = cms.bool(False),
    objects = cms.PSet(
        # -- inputs and PF --
        RawTK  = cms.VInputTag('l1tPFTracksFromL1Tracks',),
        # outputs
    ),
    copyUInts = cms.VInputTag(),
    copyFloats = cms.VInputTag(),
    copyVecUInts = cms.VInputTag(),
)
process.extraPFStuff.add(process.l1tPFTracksFromL1Tracks)

process.l1pfcandTable = cms.EDProducer("L1PFCandTableProducer",
    commonSel = cms.string("pt > 0.0 && abs(eta) < 10.0"),
    cands = cms.PSet(
    ),
    moreVariables = cms.PSet(
        # puppiWeight = cms.string("puppiWeight"),   # commented out as not a property of gen jets so raises error
        pdgId = cms.string("pdgId"),
        charge = cms.string("charge"),
    ),
)

""" Add the candidate table to the process to store candidates """
# process.l1pfgenTable = cms.EDProducer("L1PFCandTableProducer",
#     commonSel = cms.string("pt > 0.0 && abs(eta) < 10.0"),
#     cands = cms.PSet(
#         Gen = cms.InputTag("genParticles"),
#     ),
#     moreVariables = cms.PSet(
#         # puppiWeight = cms.string("puppiWeight"),   # commented out as not a property of gen jets so raises error
#         pdgId = cms.string("pdgId"),
#         charge = cms.string("charge"),
#         statusFlags = cms.string("statusFlags"),
#     ),
# )

# process.load("PhysicsTools.NanoAOD.genparticles_cff")
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer
process.l1pfgenTable = simpleCandidateFlatTableProducer.clone(
    src = cms.InputTag("genParticles"),
    name = cms.string("GenParticles"),
    doc = cms.string("gen particles"),
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

process.l1pfjetTable = cms.EDProducer("L1PFJetTableProducer",
    gen = cms.InputTag(genJets),
    commonSel = cms.string("pt > 5 && abs(eta) < 5.0"),
    drMax = cms.double(0.4 if wideJets == True else 0.2),
    minRecoPtOverGenPt = cms.double(0.1),
    jets = cms.PSet(
        Gen = cms.InputTag(genJets),
        Gen_sel = cms.string(f"pt > {str(ptCut)}"),   # str{ptCut}
    ),
    moreVariables = cms.PSet(
        nDau = cms.string("numberOfDaughters()"),
    ), 
)


"""" ADD JETS TO THE JET TABLE """
def addJets(label, tag, task, *sels):
    process.extraPFStuff.add(task)
    setattr(process.l1pfjetTable.jets, label, tag)
    for sel in sels:
        setattr(process.l1pfjetTable.jets, label+"_sel", sel)

""" SAVE JET CONSTITUENTS UP TO 128"""
def addJetConstituents(N=128):
    for i in range(N): # save a max of N daughters (unfortunately 2D arrays are not yet supported in the NanoAOD output module)
        for var in "pt", "eta", "phi", "mass", "pdgId":
            setattr(
                process.l1pfjetTable.moreVariables,    # object
                "dau%d_%s" % (i, var),    # attribute example dau0_pt
                cms.string( "? numberOfDaughters() > %d ? daughter(%d).%s : -1"  % (i, i, var) )    # value, example 1st iter:  "? numberOfDaughters() > 0 ? daughter(0).pt : -1"
                )                                                                                   # failing because number of daughters is not greater than 0 for histojets - num of daughters() returning 0
        setattr(process.l1pfjetTable.moreVariables, "dau%d_%s" % (i,"vz"), cms.string("? numberOfDaughters() > %d ? daughter(%d).%s : -1"  % (i,i,"vertex.Z")))    # Not relevant for finding seeds

""" SAVE CANDIDATES """
def saveCands(label, tag):
    setattr (process.l1pfcandTable.cands, label, cms.InputTag(tag))

############################################################################################
if wideJets == False:
    pass

if wideJets == True:
    """ GENERATOR JETS """
    gen = Jets(label="ak8Gen", tag=cms.InputTag("ak8GenJetsNoNu"), task=process.ak8GenJetsNoNuTask)
    genAK4 = Jets(label="ak4Gen", tag=cms.InputTag("ak4GenJetsNoNu"), task=process.ak4GenJetsNoNuTask)

    """ Anti-kT"""
    ak8 = Jets(label="ak8Puppi", tag=cms.InputTag("ak8PuppiJets"), task=process.ak8PuppiJetsTask)
    ak4 = Jets(label="ak4Puppi", tag=cms.InputTag("ak4PuppiJets"), task=process.ak4PuppiJetsTask)

    """ WIDE HISTOJETS (ALL TRIMMED) """
    histo18x18 = Jets(label="wideHistoPuppiEmu", tag=cms.InputTag("l1tPhase1WideJetProducer18x18", "Uncalibratedl1tPhase1WideJet18x18FromPfCandidates"), task=process.L1TPFWideJetsPhase1Task_18x18)
    histoDoubleBinSize = Jets(label="wideHistoPuppiEmuDoubleBinSize", tag=cms.InputTag("l1tPhase1WideJetProducer9x9", "Uncalibratedl1tPhase1WideJet9x9FromPfCandidates"), task=process.L1TPFWideJetsPhase1Task_9x9)

    """ SEEDED CONE 8"""
    sc8Sim = Jets(label="sc8PuppiSim", tag=cms.InputTag("l1tSC8PFL1Puppi"), task=process.L1TPFJetsTask)    # Should have jet mass
    sc8Emu = Jets(label="sc8PuppiEmu", tag=cms.InputTag("l1tSC8PFL1PuppiEmulator"), task=process.L1TPFJetsEmulationTask)    # wont have jet mass
    """ SEEDED CONE 4 """
    sc4Sim = Jets(label="sc4PuppiSim", tag=cms.InputTag("l1tSC4PFL1Puppi"), task=process.L1TPFJetsTask)
    sc4Emu = Jets(label="sc4PuppiEmu", tag=cms.InputTag("l1tSCPFL1PuppiEmulator"), task=process.L1TPFJetsEmulationTask)

    """ HISTO-SEEDED CONE 8 """
    hsc8Emu = Jets(label="hsc8PuppiEmu", tag=cms.InputTag("l1tHSC8PFL1PuppiEmu"), task=process.L1TPFHSC8JetsEmuTask)
    hsc8EmuTrimmed = Jets(label="hsc8PuppiEmuTrimmed", tag=cms.InputTag("l1tHSC8PFL1PuppiEmuTrimmed"), task=process.L1TPFHSC8JetsEmuTaskTrimmed)

    """ HISTO-SEEDED CONE 8 DOUBLE BIN SIZE """
    # Double bin size, trimmed, 9x9 mask, 1 GeV seed threshold, 0.8 cone size, deregionizer cands
    hsc8EmuDoubleBinSize = Jets(label="hsc8PuppiEmuDoubleBinSize", tag=cms.InputTag("l1tHSC8PFL1PuppiEmuDoubleBinSize"), task=process.L1TPFHSC8JetsEmuTaskDoubleBinSize)    # trimmed

addJets(*gen)
addJets(*genAK4)

addJets(*sc4Sim)
# addJets("SC8Mass120Cut", sc8Sim.tag, sc8Sim.task, cms.string("mass > 120"))
# addJets("SC8Mass110Cut", sc8Sim.tag, sc8Sim.task, cms.string("mass > 110"))
# addJets("SC8Mass100Cut", sc8Sim.tag, sc8Sim.task, cms.string("mass > 100"))
# addJets("SC8Mass90Cut",  sc8Sim.tag, sc8Sim.task, cms.string("mass > 90"))
# addJets("SC8Mass80Cut",  sc8Sim.tag, sc8Sim.task, cms.string("mass > 80"))
# addJets("SC8Mass70Cut",  sc8Sim.tag, sc8Sim.task, cms.string("mass > 70"))
# addJets("SC8Mass60Cut",  sc8Sim.tag, sc8Sim.task, cms.string("mass > 60"))
# addJets("SC8Mass50Cut",  sc8Sim.tag, sc8Sim.task, cms.string("mass > 50"))
# addJets("SC8Mass40Cut",  sc8Sim.tag, sc8Sim.task, cms.string("mass > 40"))
# addJets("SC8Mass30Cut",  sc8Sim.tag, sc8Sim.task, cms.string("mass > 30"))
# addJets("SC8Mass20Cut",  sc8Sim.tag, sc8Sim.task, cms.string("mass > 20"))
# addJets("SC8Mass10Cut",  sc8Sim.tag, sc8Sim.task, cms.string("mass > 10"))
# addJets("SC8MassNoCut",  sc8Sim.tag, sc8Sim.task                         )

# addJetConstituents(N=1)  # 128 by default as max regioniser output
saveCands("PUPPI", "l1tLayer2Deregionizer:Puppi")
# saveCands(label="GenParticles", tag = "genParticles")    # Include genParticles by default
#############################################################################################

process.p = cms.Path(
        process.ntuple + #process.content +
        process.l1pfjetTable +
        process.l1pfcandTable +
        process.l1pfgenTable
        )
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

if False:
    """ NOTE: useExternalSeeds FLAG DOESN'T WORK WITH SIM JETS, HW FLAG MUST BE TRUE """
    def wideJetStudies():
        addJets(*ak8)
        addJets(*sc8Emu)
        addJets(*hsc8Emu)

    def jetMassStudies():
        addJets("SC8Mass30Cut", sc8Sim.tag, sc8Sim.task, cms.string("mass > 30"))
        addJets("SC8Mass10Cut", sc8Sim.tag, sc8Sim.task, cms.string("mass > 10"))
        addJets("SC8MassNoCut", sc8Sim.tag, sc8Sim.task                         )
        
        addJets("AK8Mass30Cut", ak8.tag,    ak8.task,    cms.string("mass > 30"))
        addJets("AK8Mass10Cut", ak8.tag,    ak8.task,    cms.string("mass > 10"))
        addJets("AK8MassNoCut", ak8.tag,    ak8.task                            )

        addJets(*gen)

        addJetConstituents(N=1)  # 128 by default as max regioniser output
        saveCands("PUPPI", "l1tLayer2Deregionizer:Puppi")
        saveCands(label="GenParticles", tag = "genParticles")    # Include genParticles by default



