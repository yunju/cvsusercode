import FWCore.ParameterSet.Config as cms

photonDataFilter = cms.EDFilter('PhotonDataFilter',
photonLabel =cms.InputTag("selectedPatPhotons"),
jetLabel = cms.InputTag("selectedPatJetsPFlow")  #selectedPatJetsAK5PF")
 )

