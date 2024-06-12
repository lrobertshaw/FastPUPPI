import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar

import sys
from collections import namedtuple
Jets = namedtuple("Jets", "label tag task")

inputFile = str(sys.argv[-1])    # root://eoscms.cern.ch//eos/cms//store/cmst3/group/l1tr/gpetrucc/12_5_X/NewInputs125X/150223/VH_PtHat125_PU200/inputs125X_VH_PtHat125_PU200_job6.root
wideJets = True
nEvents = -1
print(f"\nRunning over file: {inputFile}\nWide jets: {wideJets}\nNumber of events: {nEvents}\n")

# Handle wide and regular jets
if wideJets == True:
    genJets = "ak8GenJetsNoNu"
    ptCut = 30
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
    inputCommands = cms.untracked.vstring("keep *", 
            "drop l1tPFClusters_*_*_*",
            "drop l1tPFTracks_*_*_*",
            "drop l1tPFCandidates_*_*_*",
            "drop l1tTkPrimaryVertexs_*_*_*")
)

process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi') # needed for HGCAL_noise_fC
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.Configuration.GenMETParticles_cff')
process.load('RecoMET.METProducers.genMetTrue_cfi')


# from RecoMET.METProducers.pfMet_cfi import pfMet

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '125X_mcRun4_realistic_v2', '')

""" DEFINE EXTRA PF STUFF AND ADD AK8 COLLECTION IF WIDEJETS == TRUE """
process.extraPFStuff = cms.Task(
    process.l1tSAMuonsGmt,
    process.l1tGTTInputProducer,
    process.l1tVertexFinderEmulator,
    process.L1TLayer1TaskInputsTask,
    process.L1TLayer1Task,
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
    #setattr(process.l1pfjetTable.jets, "AK8", cms.InputTag(ak8GenJetsNoNu))


""" DEFINE NTUPLE PRODUCER """
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

""" Add the candidate table to the process to store candidates """
process.l1pfcandTable = cms.EDProducer(
    "L1PFCandTableProducer",
    commonSel = cms.string("pt > 0.0 && abs(eta) < 10.0"),
    cands = cms.PSet(
    ),
    moreVariables = cms.PSet(
        # puppiWeight = cms.string("puppiWeight"),   # commented out as not a property of gen jets so raises error
        pdgId = cms.string("pdgId"),
        charge = cms.string("charge")
    ),
)

""" Add the jet table to the process to store jets """
process.l1pfjetTable = cms.EDProducer(
    "L1PFJetTableProducer",
    gen = cms.InputTag(genJets),
    commonSel = cms.string("pt > 5 && abs(eta) < 5.0"),
    drMax = cms.double(0.2),    # DOUBLE THE SIZE FOR WIDE CONE ????
    minRecoPtOverGenPt = cms.double(0.1),
    jets = cms.PSet(
        Gen = cms.InputTag(genJets),
        Gen_sel = cms.string(f"pt > {str(ptCut)}"),
    ),
    moreVariables = cms.PSet(
        nDau = cms.string("numberOfDaughters()"),
    ),
)

"""" ADD JETS TO THE JET TABLE """
def addJets(label, tag, task):
    process.extraPFStuff.add(task)
    setattr(process.l1pfjetTable.jets, label, tag)

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

##############################################################################################################################################################

if wideJets == False:
    """ HISTOJETS 9X9 """
    # Untrimmed
    histo9x9Emu = Jets(label="histo9x9Emu", tag=cms.InputTag('l1tPhase1PuppiJetProducer9x9', "Uncalibratedl1tPhase1JetFromPuppiCandidates"), task=process.L1TPuppiJetsPhase1Task_9x9)
    # Trimmed
    histo9x9trimmedEmu = Jets(label="histo9x9trimmedEmu", tag=cms.InputTag('l1tPhase1PuppiJetProducer9x9trimmed', "Uncalibratedl1tPhase1TrimmedJetFromPuppiCandidates"), task=process.L1TPuppiJetsPhase1Task_9x9trimmed)

    """ HISTOJETS 7x7 """
    histo7x7Emu = Jets(label="histo7x7Emu", tag=cms.InputTag('l1tPhase1PuppiJetProducer7x7', "Uncalibratedl1tPhase1Jet7x7FromPuppiCandidates"), task=process.L1TPuppiJetsPhase1Task_7x7)
    histo7x7trimmedEmu = Jets(label="histo7x7trimmedEmu", tag=cms.InputTag('l1tPhase1PuppiJetProducer7x7trimmed', "Uncalibratedl1tPhase1Trimmed7x7JetFromPuppiCandidates"), task=process.L1TPuppiJetsPhase1Task_7x7trimmed)

    """ SEEDED CONE 4 """
    sc4Emu = Jets(label="sc4PuppiEmu", tag=cms.InputTag("l1tSCPFL1PuppiEmulator"), task=process.L1TPFJetsTask)

    """ HISTO-SEEDED CONE 4 """
    # Untrimmed
    hsc4Emu = Jets(label="hsc4PuppiEmu", tag=cms.InputTag("l1tHSCPFL1PuppiEmulator"), task=process.L1TPFHSCJetsEmuTask)
    # Trimmed
    hsc4EmuTrimmed = Jets(label="hsc4PuppiEmuTrimmed", tag=cms.InputTag("l1tHSCPFL1PuppiEmulatorTrimmed"), task=process.L1TPFHSCJetsEmuTaskTrimmed)

    """ HISTO-SEEDED CONE 4 EMULATOR JETS WITH DIFFERENT MASK SIZES """
    # Untrimmed
    hsc4Emu9x9mask = Jets(label="hsc4PuppiEmu9x9", tag=cms.InputTag("l1tHSCPFL1Puppi9x9Emu"), task=process.l1tHSCPFL1Puppi9x9EmuTask)
    hsc4Emu7x7mask = Jets(label="hsc4PuppiEmu7x7", tag=cms.InputTag("l1tHSCPFL1Puppi7x7Emu"), task=process.l1tHSCPFL1Puppi7x7EmuTask)
    hsc4Emu5x5mask = Jets(label="hsc4PuppiEmu5x5", tag=cms.InputTag("l1tHSCPFL1Puppi5x5Emu"), task=process.l1tHSCPFL1Puppi5x5EmuTask)
    hsc4Emu3x3mask = Jets(label="hsc4PuppiEmu3x3", tag=cms.InputTag("l1tHSCPFL1Puppi3x3Emu"), task=process.l1tHSCPFL1Puppi3x3EmuTask)
    # Trimmed
    hsc4Emu9x9maskTrimmed = Jets(label="hsc4PuppiEmu9x9Trimmed", tag=cms.InputTag("l1tHSCPFL1Puppi9x9EmuTrimmed"), task=process.l1tHSCPFL1Puppi9x9EmuTaskTrimmed)
    hsc4Emu7x7maskTrimmed = Jets(label="hsc4PuppiEmu7x7Trimmed", tag=cms.InputTag("l1tHSCPFL1Puppi7x7EmuTrimmed"), task=process.l1tHSCPFL1Puppi7x7EmuTaskTrimmed)
    hsc4Emu5x5maskTrimmed = Jets(label="hsc4PuppiEmu5x5Trimmed", tag=cms.InputTag("l1tHSCPFL1Puppi5x5EmuTrimmed"), task=process.l1tHSCPFL1Puppi5x5EmuTaskTrimmed)
    hsc4Emu3x3maskTrimmed = Jets(label="hsc4PuppiEmu3x3Trimmed", tag=cms.InputTag("l1tHSCPFL1Puppi3x3EmuTrimmed"), task=process.l1tHSCPFL1Puppi3x3EmuTaskTrimmed)

    """ HISTO-SEEDED CONE 4 EMULATOR JETS WITH DIFFERENT MASK AND SEED SIZES """
    # Untrimmed
    hsc4Emu9x9mask1x1seed = Jets(label="hsc4PuppiEmu9x9Mask1x1Seed", tag=cms.InputTag("l1tHSCPFL1Puppi9x9Seed1x1Emu"), task=process.l1tHSCMaskSeedSizesJetsTask)
    hsc4Emu9x9mask3x3seed = Jets(label="hsc4PuppiEmu9x9Mask3x3Seed", tag=cms.InputTag("l1tHSCPFL1Puppi9x9Seed3x3Emu"), task=process.l1tHSCMaskSeedSizesJetsTask)
    hsc4Emu7x7mask1x1seed = Jets(label="hsc4PuppiEmu7x7Mask1x1Seed", tag=cms.InputTag("l1tHSCPFL1Puppi7x7Seed1x1Emu"), task=process.l1tHSCMaskSeedSizesJetsTask)
    hsc4Emu7x7mask3x3seed = Jets(label="hsc4PuppiEmu7x7Mask3x3Seed", tag=cms.InputTag("l1tHSCPFL1Puppi7x7Seed3x3Emu"), task=process.l1tHSCMaskSeedSizesJetsTask)
    # Trimmed
    hsc4Emu9x9mask1x1seedTrimmed = Jets(label="hsc4PuppiEmu9x9Mask1x1SeedTrimmed", tag=cms.InputTag("l1tHSCPFL1Puppi9x9Seed1x1EmuTrimmed"), task=process.l1tHSCMaskSeedSizesJetsTaskTrimmed)
    hsc4Emu9x9mask3x3seedTrimmed = Jets(label="hsc4PuppiEmu9x9Mask3x3SeedTrimmed", tag=cms.InputTag("l1tHSCPFL1Puppi9x9Seed3x3EmuTrimmed"), task=process.l1tHSCMaskSeedSizesJetsTaskTrimmed)
    hsc4Emu7x7mask1x1seedTrimmed = Jets(label="hsc4PuppiEmu7x7Mask1x1SeedTrimmed", tag=cms.InputTag("l1tHSCPFL1Puppi7x7Seed1x1EmuTrimmed"), task=process.l1tHSCMaskSeedSizesJetsTaskTrimmed)
    hsc4Emu7x7mask3x3seedTrimmed = Jets(label="hsc4PuppiEmu7x7Mask3x3SeedTrimmed", tag=cms.InputTag("l1tHSCPFL1Puppi7x7Seed3x3EmuTrimmed"), task=process.l1tHSCMaskSeedSizesJetsTaskTrimmed)

if wideJets == True:
    """ GENERATOR JETS """
    #ak8 = Jets(label="AK8", tag=cms.InputTag("ak8GenJetsNoNu"), task=process.ak8GenJetsNoNuTask)

    """ SEEDED CONE 8"""
    sc8Sim = Jets(label="sc8PuppiSim", tag=cms.InputTag("l1tSCPFL1PuppiSimWide"), task=process.L1TPFWideJetsSimTask)
    sc8Emu = Jets(label="sc8PuppiEmu", tag=cms.InputTag("l1tSCPFL1PuppiEmuWide"), task=process.L1TPFWideJetsEmuTask)

    """ HISTO-SEEDED CONE 8 """
    hsc8Emu = Jets(label="hsc8PuppiEmu", tag=cms.InputTag("l1tHSCPFL1PuppiEmuWide"), task=process.L1TPFHSCWideJetsEmuTask)
    hsc8EmuTrimmed = Jets(label="hsc8PuppiEmuTrimmed", tag=cms.InputTag("l1tHSCPFL1PuppiWideEmulatorTrimmed"), task=process.L1TPFHSCWideJetsEmuTaskTrimmed)


addJets(*sc8Emu)
addJets(*hsc8Emu)
addJets(*hsc8EmuTrimmed)


# addJetConstituents(N=128)  # 128 by default

# saveCands(label="l1tLayer2DeregionizerPUPPI", tag="l1tLayer2Deregionizer:Puppi")    # save emu PUPPI cands
# saveCands(label="l1tLayer1PUPPI", tag="l1tLayer1:Puppi")    # save sim PUPPI cands
# saveCands(label="GenParticles", tag = "genParticles")    # Include genParticles by default

##############################################################################################################################################################

""" END PROCESS AND OUTPUT """
process.p = cms.Path(
    process.ntuple +
    process.l1pfjetTable +
    process.l1pfcandTable
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
    def regJetsStudies():
        # Untrimmed
        addJets(*histo9x9Emu)
        addJets(*sc4Emu)
        addJets(*hsc4Emu)
        # Trimmed
        addJets(*histo9x9trimmedEmu)
        addJets(*sc4Emu)
        addJets(*hsc4EmuTrimmed)

    def maskSizeStudies():
        # Untrimmed
        addJets(*hsc4Emu9x9mask)
        addJets(*hsc4Emu7x7mask)
        addJets(*hsc4Emu5x5mask)
        addJets(*hsc4Emu3x3mask)
        # Trimmed
        addJets(*hsc4Emu9x9maskTrimmed)
        addJets(*hsc4Emu7x7maskTrimmed)
        addJets(*hsc4Emu5x5maskTrimmed)
        addJets(*hsc4Emu3x3maskTrimmed)

    def seedSizeStudies():
        pass

    def maskAndSeedStudies():
        # Untrimmed
        addJets(*hsc4Emu9x9mask1x1seed)
        addJets(*hsc4Emu9x9mask3x3seed)
        addJets(*hsc4Emu7x7mask1x1seed)
        addJets(*hsc4Emu7x7mask3x3seed)
        # Trimmed
        addJets(*hsc4Emu9x9mask1x1seedTrimmed)
        addJets(*hsc4Emu9x9mask3x3seedTrimmed)
        addJets(*hsc4Emu7x7mask1x1seedTrimmed)
        addJets(*hsc4Emu7x7mask3x3seedTrimmed)


    def wideJetStudies():
        addJets(*sc8Emu)
        addJets(*hsc8Emu)
