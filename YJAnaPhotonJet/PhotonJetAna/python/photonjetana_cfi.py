import FWCore.ParameterSet.Config as cms

photonjetana = cms.EDAnalyzer('PhotonJetAna',
                      doGenParticles = cms.bool(True),
                      VertexProducerPY  = cms.InputTag("offlinePrimaryVertices"),                     
                      BeamSpotProducerPY =   cms.InputTag("offlineBeamSpot"),  
                      JetProducerPY =   cms.InputTag("cleanPatJets"),
                      triggerEventProducerPY =   cms.InputTag("patTriggerEvent"),  
                      genParticlesProducerPY =   cms.InputTag("genParticles"),
                      YJJetProducerPY        =   cms.InputTag("cleanPatJets"),
                      PhotonProducerPY       =   cms.InputTag("selectedPatPhotons"),  
                      ebReducedRecHitCollectionPY =  cms.InputTag("reducedEcalRecHitsEB"),
                      eeReducedRecHitCollectionPY =  cms.InputTag("reducedEcalRecHitsEE"),                      
                      
)
