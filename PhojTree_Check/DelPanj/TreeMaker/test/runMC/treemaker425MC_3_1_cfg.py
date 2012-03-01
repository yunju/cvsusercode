import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
 ## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

# for 2011 42  data
#process.GlobalTag.globaltag = cms.string('GR_R_42_V12::All')
# for 2011 42  mc
process.GlobalTag.globaltag = cms.string('START42_V17::All')



process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')


#pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
## need for run data

#runOnData(process, ['All'], outputInProcess = False)



# add PFMet
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

process.source = cms.Source("PoolSource",
   # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
     #rfio:/castor/cern.ch/user/l/lovedeep/Spring10/Prod2_WJets-madgraph/patLayer1_fromAOD_PF2PAT_Prod2_90_1_XSl.root"
     #$inputFileNames
     #'file:/hdfs/store/user/kaur/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauolaSummer11-PU_S4_START42_V11-v1_DelPanj_Jul29-Electrons_PF2PAT_onMC/Electrons_PF2PAT_onMC-002FE237-B09C-E011-B7B1-0022199305B1.root'
      #'file:/castor/cern.ch/user/a/anil/patTupleZee_PF2PAT.root'
      #'file:/scratch/yunju/PhoJAna428_2011A/CMSSW_4_2_8/src/48FA4B21-DA7F-E011-8903-003048CFB40C.root'
#      'file:/data4/yunju/Validation/QCD_Pt-120to170_TuneZ2_7TeV_pythia6-AODSIM-PU_S4_START42_V11-v.root'
    'file:/data4/yunju/NewtreeT/VaAod/pho170val.root' 
  )
)
#process.load("DelPanj.TreeMaker.patTuples_ZmmMadgraph_cfi")
baseJetSel = cms.PSet(
  #Jets=cms.InputTag("selectedPatJetsPFlow")
   Jets=cms.InputTag("selectedPatJetsAK5PF")
)

process.demo = cms.EDAnalyzer('TreeMaker',
   
   fillEventInfo_ = cms.bool(True),
   fillGenInfo_   = cms.bool(True),
   fillMuonInfo_  = cms.bool(False),
   fillElecInfo_  = cms.bool(False),
   fillJetInfo_   = cms.bool(True),
   fillMetInfo_   = cms.bool(False),
   fillTrigInfo_  = cms.bool(True),
   fillPhotInfo_  = cms.bool(True),	
   genPartLabel=cms.InputTag("genParticles"),
   patMuons=cms.InputTag("selectedPatMuons"),
   patElectrons = cms.InputTag("selectedPatElectronsPFlow"),
   Jets=cms.InputTag("selectedPatJetsAK5PF"),
   patMet=cms.InputTag("patMETs"), 
   beamSpotLabel=cms.InputTag("offlineBeamSpot"),
   VertexProducerPY  = cms.InputTag("offlinePrimaryVertices"),
   BeamSpotProducerPY =   cms.InputTag("offlineBeamSpot"),
   rho25PY = cms.InputTag("kt6PFJets25:rho"),
   rho44PY = cms.InputTag("kt6PFJets44:rho"),
   photonLabel =cms.InputTag("selectedPatPhotons"),
   reducedEcalRecHitsEBLabel =  cms.InputTag("reducedEcalRecHitsEB"),
   reducedEcalRecHitsEELabel =  cms.InputTag("reducedEcalRecHitsEE"),
   genParticlesProducerPY =cms.InputTag("genParticles") ,
   isMCPY = cms.bool(True),
   HltKeyWordPY =cms.string("HLT_Photon"),   
   usePFlow = cms.bool(False),
   patJetPfAk05 = baseJetSel,
   #outFileName=cms.string('$outputFileName')
   outFileName=cms.string('Phojtree.root')
)

from DelPanj.TreeMaker.zeeFilter_cff import *
process.zeefil = zeeFilter
process.zeefil.leadElecPset_.ptx= cms.double(20)
process.zeefil.subLeadElecPset_.ptx= cms.double(20)

from DelPanj.TreeMaker.photonFilter_cff import *
process.photonfil = photonFilter
process.photonfil.GammaPtMinPY=cms.double(40)

#process.p = cms.Path(process.zeefil*process.demo)
process.p = cms.Path(
                     process.fjSequence25*
                     process.fjSequence44*
                     process.patDefaultSequence *
                     photonFilter*  
                     process.demo 
                     )


#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
 



