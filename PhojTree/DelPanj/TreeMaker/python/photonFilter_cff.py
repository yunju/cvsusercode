import FWCore.ParameterSet.Config as cms

photonFilter = cms.EDFilter('PhotonFilter',
photonLabel =cms.InputTag("selectedPatPhotons"),
jetLabel = cms.InputTag("selectedPatJetsAK5PF")
 )

