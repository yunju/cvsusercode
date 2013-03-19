import FWCore.ParameterSet.Config as cms

process = cms.Process("COMBINE")
#runOnMC = False # to run on DATA
runOnMC = True

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
#  'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/VBF_HToZZTo2L2Q_M-200_8TeV-powheg-pythia6_h2l2qSkimData_1_1_ClF.root'
#   'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/GluGluToHToZZTo2L2Q_M-200_8TeV-powheg-pythia6_SkimPAT_H2l2q_523_v3_l_h2l2qSkimData_1_1_JVm.root'
#'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/DYJetsToLL_M-10To50_h2l2qSkimData_51_1_2MJ.root',
#  'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/DYJetsToLL_M-10To50_h2l2qSkimData_205_1_stI.root',
#'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/ZZ_TuneZ2star_8TeV_pythia6_tauolah2l2qSkimData_382_1_foR.root' 
#'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/DYJetsToLL_M-10To50_h2l2qSkimData_104_1_sZg.root',
#'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/DYJetsToLL_M-10To50filter_8TeV-madgraph_h2l2qSkimData_3_1_v1i.root',
#'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/DYJetsToLL_M-10To50_h2l2qSkimData_15_1_Jz1.root',
#'file:/scratch/yunju/2kTW/CMSSW_5_2_5/src/DYJetsToLL_M-10To50_h2l2qSkimData_219_1_RZs.root'
#'root://eoscms//eos/cms/store/mc/Summer12/GluGluToHToZZTo2L2Q_M-200_8TeV-powheg-pythia6/AODSIM/PU_S7_START1B1A-3AB8-E111-8D00-003048FFCB8C.root'
#'file:/data4/yunju/VBF2012/skim2l2q/SkimPAT_H2l2q_523_v3_l_GluGluToHToZZTo2L2Q_M-200h2l2qSkimData_22_1_2KD.root'
#'file:/home/yunju/cdxfe_Higgs/Test_525skim/skimsample/YJTest_h2l2qSkimData.root'
'file:/data4/yunju/VBF2012/SkimPAT_H2l2q_533_v1_GluGluToHToZZTo2L2Q_M-200_h2l2qSkimData_2_1_5SE_2nd.root'
   )
 )

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1))

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB062012")
process.load("RecoBTag.PerformanceDB.BTagPerformanceDB062012")


process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#baseJetSel = cms.PSet(
#  Jets=cms.InputTag("cleanPatJetsNoPUIsoLept")
#   JetsPY=cms.InputTag("customPFJets")
#)

#baseJetSel = cms.PSet(
#  Jets=cms.InputTag("pfJetsAK5")
#)


from DelPanj.TreeMaker.eSel2012HZZ_cff import *
from DelPanj.TreeMaker.muSelYJVBF_cff import *


process.tree1 = cms.EDAnalyzer(
	'TreeMaker',
	fillPUweightInfo_ = cms.bool(False),
	fillEventInfo_ = cms.bool(False),
	fillGenInfo_   = cms.bool(False),
	fillMuonInfo_  = cms.bool(True),
	fillElecInfo_  = cms.bool(True),
	fillElecIsoInfo_ = cms.bool(False),
	fillJetInfo_   = cms.bool(True),
	fillMetInfo_   = cms.bool(False),
	fillTrigInfo_  = cms.bool(False),
	fillPhotInfo_  = cms.bool(False),
	fillZJetPlant_ = cms.bool(False),
	fillZZInfo_    = cms.bool(False),
        
        eleRhoIso = cms.InputTag("kt6PFJetsForIso","rho"),
        muoRhoIso = cms.InputTag("kt6PFJetsCentralNeutral", "rho"),
        
        e2012IDSet  =  eSel2012HZZ,
        mu2012IDSet = muSelYJVBF,
        
        hzzeejjTag = cms.InputTag("hzzeejj:h"),
        hzzmmjjTag = cms.InputTag("hzzmmjj:h"),
	
        genPartLabel=cms.InputTag("genParticles"),
	
        patMuonsPY=cms.InputTag("userDataSelectedMuons"),
	patElectronsPY = cms.InputTag("userDataSelectedElectrons"),
        JetsPY=cms.InputTag("customPFJetsNoPUSub"), 
        

        JetForElecTreePY =cms.InputTag("kt6PFJets"),
        leadElecPset_ = eSel2012HZZ,
	patMet=cms.InputTag("patMETs"),
	beamSpotLabel=cms.InputTag("offlineBeamSpot"),
	outFileName=cms.string('outputFileName.root'),
	
)



process.counter_original = cms.EDAnalyzer('YJEventCounter',
   instance = cms.int32(1) 
)


process.TFileService = cms.Service("TFileService",
      #fileName = cms.string("DYJ200YJTest_v2.root"),
    #  fileName = cms.string("ZZ200YJTest_v2.root"),
      #fileName = cms.string("GGH200YJTest_v2.root"),
      #fileName = cms.string("Mjj_VBF200YJTest_v2.root"),
      fileName = cms.string("VBFTree_MC_53X.root"),




       closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(
	process.counter_original*
	process.tree1##Trigger Applied.
	)
  
 



