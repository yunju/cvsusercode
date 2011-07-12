import FWCore.ParameterSet.Config as cms

process = cms.Process("PhotonJetAna")

process.load("FWCore.MessageService.MessageLogger_cfi")
## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# for data
process.GlobalTag.globaltag = cms.string('START41_V0::All')

process.load("Configuration.StandardSequences.MagneticField_cff")

process.TFileService = cms.Service("TFileService", fileName = cms.string('Photonjetana.root'))


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.coreTools import *

# Jets
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJets14 = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets14.Rho_EtaMax = cms.double(1.4442)
process.fjSequence14 = cms.Sequence( process.kt6PFJets14 )

process.kt6PFJets25 = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)
process.fjSequence25 = cms.Sequence( process.kt6PFJets25 )

process.kt6PFJets30 = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets30.Rho_EtaMax = cms.double(3.0)
process.fjSequence30 = cms.Sequence( process.kt6PFJets30 )

process.kt6PFJets44 = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
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
      'file:/data2/yunju/Summer2011/PhoJet_15AODSIM.root'

    )
)

 # load the coreTools of PAT

#from PhysicsTools.PatAlgos.tools.metTools import *
#addTcMET(process, 'TC')
#addPfMET(process, 'PF')



process.load("YJAnaPhotonJet.PhotonJetAna.photonjetana_cfi") 


process.photonjetana.YJJetProducerPY = cms.InputTag("selectedPatJetsAK5PF")


process.p = cms.Path(
                      process.fjSequence14*
                      process.fjSequence25*
                      process.fjSequence30*
                      process.fjSequence44*
                      process.patDefaultSequence *
                      process.photonjetana)
