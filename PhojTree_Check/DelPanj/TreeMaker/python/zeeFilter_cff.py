import FWCore.ParameterSet.Config as cms
from DelPanj.TreeMaker.eSelLvdp2011_cff import eSelLvdp2011

zeeFilter = cms.EDFilter('ZeeFilter',
  patElectrons = cms.InputTag("selectedPatElectronsPFlow"), 
  leadElecPset_ = eSelLvdp2011,
  subLeadElecPset_ = eSelLvdp2011 
)

