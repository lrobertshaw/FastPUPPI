import FWCore.ParameterSet.Config as cms

from FastPUPPI.NtupleProducer.runPerformanceNTuple import *

noResp();
addGenLep([11,22]);
addTkEG(); # done
addHGCalTPs() # done
addStaEG() # done
noResp()
addDecodedTk(); # done
addEGCrystalClusters(); # done
# addSeededConeJets(); # done
addTkJets();
addStaMu();
addPFLep();


process.maxEvents.input = cms.untracked.int32(100)
process.source.fileNames = cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/cmst3/group/l1tr/FastPUPPI/14_0_X/fpinputs_131X/v9a/DoubleElectron_FlatPt-1To100_PU200/inputs131X_1.root',
'root://eoscms.cern.ch//eos/cms/store/cmst3/group/l1tr/FastPUPPI/14_0_X/fpinputs_131X/v9a/DoubleElectron_FlatPt-1To100_PU200/inputs131X_10.root')
process.outnano.fileName = cms.untracked.string('perfNano.root')

