import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
 ## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

#pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
## need for run data

runOnData(process, ['All'], outputInProcess = False)

#vertexing
from PhysicsTools.PatAlgos.patTemplate_cfg import *
del process.source

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

# for 2011 42 data
process.GlobalTag.globaltag = cms.string('GR_R_42_V23::All')

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.patJetCorrFactors.rho = cms.InputTag('kt6PFJetsPFlow','rho')

process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)
process.fjSequence25 = cms.Sequence( process.kt6PFJets25 )

##process.kt6PFJets44 = process.kt6PFJets.clone( doRhoFastjet = True )
##process.kt6PFJets44.Rho_EtaMax = cms.double(4.4)
##process.kt6PFJets44.Ghost_EtaMax = cms.double(5.0)
##process.fjSequence44 = cms.Sequence( process.kt6PFJets44 )

#pf2pat
from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PFlow"
usePF2PAT(process,runPF2PAT=True,jetAlgo='AK5',runOnMC=True,postfix = postfix,jetCorrections=('AK5PFchs',['L1FastJet','L2Relative','L3Absol\
ute','L2L3Residual']))
process.pfPileUpPFlow.Enable = True
process.pfPileUpPFlow.Vertices = 'goodOfflinePrimaryVertices'
process.pfJetsPFlow.doAreaFastjet = True
process.pfJetsPFlow.doRhoFastjet =  False
process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)

process.patJetCorrFactorsPFlow.rho = cms.InputTag("kt6PFJetsPFlow","rho")

# Compute the mean pt per unit area (rho) from the PF chs inputs
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsPFlow = kt4PFJets.clone(
    rParam = cms.double(0.6),
    src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
   # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/data/Run2011A/Photon/AOD/May10ReReco-v1/0000/00372EF6-C37B-E011-983B-001A64789E0C.root'
   )
)
#process.load("DelPanj.TreeMaker.patTuples_ZmmMadgraph_cfi")
baseJetSel = cms.PSet(
    Jets=cms.InputTag("selectedPatJetsPFlow")
    #Jets=cms.InputTag("selectedPatJetsAK5PF")
)

process.demo = cms.EDAnalyzer('TreeMaker',
   
   fillEventInfo_ = cms.bool(True),
   fillGenInfo_   = cms.bool(False),
   fillMuonInfo_  = cms.bool(False),
   fillElecInfo_  = cms.bool(False),
   fillJetInfo_   = cms.bool(True),
   fillMetInfo_   = cms.bool(False),
   fillTrigInfo_  = cms.bool(True),
   fillPhotInfo_  = cms.bool(True),	
   genPartLabel=cms.InputTag("genParticles"),
   patMuons=cms.InputTag("selectedPatMuons"),
   patElectrons = cms.InputTag("selectedPatElectronsPFlow"),
   Jets=cms.InputTag("selectedPatJetsPFlow"), # selectedPatJetsAK5PF"),
   patMet=cms.InputTag("patMETs"), 
   beamSpotLabel=cms.InputTag("offlineBeamSpot"),
   VertexProducerPY  = cms.InputTag("offlinePrimaryVertices"),
   BeamSpotProducerPY =   cms.InputTag("offlineBeamSpot"),
   rho25PY = cms.InputTag("kt6PFJets25:rho"),
   rho44PY = cms.InputTag("kt6PFJetsPFlow:rho"), #44:rho"),
   photonLabel =cms.InputTag("selectedPatPhotons"),
   reducedEcalRecHitsEBLabel =  cms.InputTag("reducedEcalRecHitsEB"),
   reducedEcalRecHitsEELabel =  cms.InputTag("reducedEcalRecHitsEE"),
   genParticlesProducerPY =cms.InputTag("genParticles") ,
   isMCPY = cms.bool(False),
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
process.photonfil.GammaPtMinPY =cms.double(40)

from DelPanj.TreeMaker.photonDataFilter_cff import *
process.photondatafil = photonDataFilter
process.photondatafil.GammaPtMinPY =cms.double(40)
process.photondatafil.GammaEtaMaxPY =cms.double(2.6)

# Add the KT6 producer to the sequence
getattr(process,"patPF2PATSequence"+postfix).replace(
    getattr(process,"pfNoElectron"+postfix), getattr(process,"pfNoElectron"+postfix) * process.kt6PFJets25 * process.kt6PFJetsPFlow
    )
process.patseq = cms.Sequence(
    process.goodOfflinePrimaryVertices*
    getattr(process,"patPF2PATSequence"+postfix)
    )

#Can add to the path to see what is included in the event and what's not --- used for debugging purposes
process.content = cms.EDAnalyzer("EventContentAnalyzer")

#beam scraping filter
process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(True),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )

#gets rid out the patTuple.root output normally given when running PAT
del process.outpath

#process.p = cms.Path(process.zeefil*process.demo)
process.p = cms.Path(
                     process.noscraping*
#                     process.primaryVertexFilter*
#                     process.fjSequence25*
#                     process.fjSequence44*
                     process.patseq *
                     process.patDefaultSequence *
                     process.photondatafil*  
                     process.demo 
                     )


#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
 



