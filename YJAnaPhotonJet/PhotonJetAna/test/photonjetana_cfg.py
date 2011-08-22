import FWCore.ParameterSet.Config as cms

process = cms.Process("PhotonJetAna")

process.load("FWCore.MessageService.MessageLogger_cfi")
## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

# for 2011 42  data
process.GlobalTag.globaltag = cms.string('GR_R_42_V12::All')

process.load("Configuration.StandardSequences.MagneticField_cff")

process.TFileService = cms.Service("TFileService", fileName = cms.string('YJPhotonjetana2011.root'))


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
## need for run data
runOnData(process, ['All'], outputInProcess = False)

# Jets
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

#process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.load('RecoJets.Configuration.RecoPFJets_cff')

process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)
process.fjSequence25 = cms.Sequence( process.kt6PFJets25 )

process.kt6PFJets44 = process.kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets44.Rho_EtaMax = cms.double(4.4)
process.kt6PFJets44.Ghost_EtaMax = cms.double(5.0)
process.fjSequence44 = cms.Sequence( process.kt6PFJets44 )



process.patJetCorrFactors.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets44','rho')





# add PFMet
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')



from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,
                  cms.InputTag('ak5PFJets'),
                  'AK5', 'PF',
                  doJTA        = True,
                  doBTagging   = True,
                  jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
                  doType1MET   = False,
                  doL1Cleaning = False,
                  doL1Counters = False,
                  genJetCollection = cms.InputTag("ak5GenJets"),
                  doJetID      = False
                  )




process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
      'file:/scratch/yunju/PhoJAna428_2011A/CMSSW_4_2_8/src/48FA4B21-DA7F-E011-8903-003048CFB40C.root'

    )
)

 # load the coreTools of PAT

#from PhysicsTools.PatAlgos.tools.metTools import *
#addTcMET(process, 'TC')
#addPfMET(process, 'PF')



process.load("YJAnaPhotonJet.PhotonJetAna.photonjetana_cfi") 


process.photonjetana.YJJetProducerPY = cms.InputTag("selectedPatJetsAK5PF")


process.p = cms.Path(
                      #process.fjSequence14*
                      process.fjSequence25*
                     # process.fjSequence30*
                      process.fjSequence44*
                      process.patDefaultSequence *
                      process.photonjetana)
