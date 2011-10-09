#ifndef __PHOTON_TREE_H_
#define __PHOTON_TREE_H_

/*
Log:
Sep 10, 2011
Anil Singh: Empty template created. 
*/

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "DelPanj/TreeMaker/interface/baseTree.h"

 #include "DataFormats/EcalDetId/interface/EBDetId.h"
  #include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
 #include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

 #include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
 #include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
 #include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
 #include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
 #include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
 #include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
 #include "DataFormats/ParticleFlowReco/interface/PFBlockElementSuperCluster.h"
// #include "DataFormats/EgammaCandidates/interface/Photon.h"
//MPA

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
 
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"

//Trigger DataFormats
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"

#include "DataFormats/Common/interface/TriggerResults.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CommonTools/Utils/interface/PtComparator.h"
//#include "RecoEgamma/EgammaTools/interface/ConversionLikelihoodCalculator.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include <Math/VectorUtil.h>
#include <TLorentzVector.h>


using namespace std;
using namespace edm;

class photonTree : public baseTree{

 public:
  photonTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~photonTree();
  void Fill(const edm::Event& iEvent,const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();

 private:
  photonTree(){};
  bool usePFObjects_;
  edm::InputTag photonLabel_;
  edm::InputTag reducedEcalRecHitsEBLabel_;
  edm::InputTag reducedEcalRecHitsEELabel_;
  edm::InputTag genParticlesProducerC;
  edm::InputTag rho25C;
  edm::InputTag rho44C;

  int          pdgId_; 
  std::vector<int> otherPdgIds_; 
  int muInGenCaloIso_;
  
  
  virtual Float_t getGenCalIso(edm::Handle<reco::GenParticleCollection> handle,
                              reco::GenParticleCollection::const_iterator thisPho, const Float_t dRMax=0.4, bool removeMu=true, bool removeNu=false);

  virtual Float_t getGenTrkIso(edm::Handle<reco::GenParticleCollection> handle,
                              reco::GenParticleCollection::const_iterator thisPho, const Float_t dRMax=0.4);


  bool isMC_;
  //variables which would become branches
  double rho25_;
  double rho44_;
  double nPhoton_;
  std::vector<double> photonPt_;
  std::vector<double> photonEta_;
  std::vector<double> photonPhi_;
  std::vector<double> photonEt_;
  std::vector<double> photonEnergy_;
  std::vector<double> photonPx_;
  std::vector<double> photonPy_;
  std::vector<double> photonPz_;
  std::vector<double> photonR9_;
  std::vector<double> photonPhiWidth_;
  std::vector<double> photonEtaWidth_;
  std::vector<double> photonScPhi_;
  std::vector<double> photonScEta_;
  std::vector<double> photonESRatio_;
  std::vector<double> photonSigmaIetaIeta_;
  std::vector<double> photonSigmaIphiIphi_;
  std::vector<double> photonSeedTime_;
  std::vector<double> photonseedSeverity_;
  std::vector<double> photonE2overe9_;
  std::vector<double> photonhadronicOverEm_;
  std::vector<double> photonecalRecHitSumEtConeDR04_;
  std::vector<double> photonhcalTowerSumEtConeDR04_;
  std::vector<double> photontrkSumPtHollowConeDR04_;
  std::vector<double> photonhasPixelSeed_;
  double  MCpthat_;
  std::vector<bool>  PhoisGenMatched_;
  std::vector<double>   PhogenMomId_; 
  std::vector<double>   PhogenGrandMomId_;
  std::vector<double>   PhogenNSiblings_;
  std::vector<double>   PhogenMatchedE_;
  std::vector<double>   PhogenMatchedPx_;
  std::vector<double>   PhogenMatchedPy_;
  std::vector<double>   PhogenMatchedPz_;
  std::vector<double>   PhogenMatchedPt_;
  std::vector<double>   PhogenMatchedEta_;
  std::vector<double>   PhogenMatchedPhi_;
  std::vector<double>   PhogenCalIsoDR04_;
  std::vector<double>   PhogenTrkIsoDR04_;
  std::vector<double>   PhogenIsoDR04_;


 

 
};

#endif

