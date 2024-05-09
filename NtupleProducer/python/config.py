import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar

import sys
from collections import namedtuple
Jets = namedtuple("Jets", "label tag task")

inputFile = sys.argv[-1]    # root://eoscms.cern.ch//eos/cms//store/cmst3/group/l1tr/gpetrucc/12_5_X/NewInputs125X/150223/VH_PtHat125_PU200/inputs125X_VH_PtHat125_PU200_job6.root
print(inputFile)

nEvents = -1    # -1 is all events
wideJets = True

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
    setattr(process.l1pfjetTable.jets, label, cms.InputTag(tag))

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
AK8 = Jets(label="AK8", tag="ak8GenJetsNoNu", task=process.ak8GenJetsNoNuTask)
SC8 = Jets(label="SC8", tag="l1tSCPFL1PuppiEmulator", task=process.L1TPFJetsTask)
wideHSC = Jets(label="wideHSC", tag="l1t17x17HistoSeedsSCPFL1PuppiEmulator", task=process.L1TPFHistoSeedJetsTask)
wideHSC3x3 = Jets(label="wideHSC3x3", tag="l1t17x17Histo3x3SeedsSCPFL1PuppiEmulator", task=process.L1TPFHisto3x3SeedJetsTask)
trimmedWideHSC = Jets(label="trimmedWideHSC", tag="l1t17x17TrimmedHistoSeedsSCPFL1PuppiEmulator", task=process.L1TPFTrimmedHistoSeedJetsTask)
wide15x15HSC3x3 = Jets(label="wide15x15HSC3x3", tag="l1t15x15Histo3x3SeedsSCPFL1PuppiEmulator", task=process.L1TPF15x15Histo3x3SeedJetsTask)

addJets(*AK8)
addJets(*SC8)
addJets(*wideHSC)
# addJets(*wideHSC3x3)
# addJets(*trimmedWideHSC)
# addJets(*wide15x15HSC3x3)

addJetConstituents(N=128)  # 128 by default

saveCands(label="l1tLayer2Deregionizer", tag="l1tLayer2Deregionizer:Puppi")    # save PUPPI cands
saveCands(label="GenParticles", tag = "genParticles")    # Include genParticles by default
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