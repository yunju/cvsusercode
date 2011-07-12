#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//from gg

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

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


 #include <memory>
 #include <fstream>
 #include <map>


 #include "TTree.h"
 #include "TH1F.h"
 #include "TH2F.h"
 #include "TString.h"





using namespace edm;
using namespace std;
//using namespace reco;
using namespace pat;
using namespace pat::helper;

const Int_t maxP = 500;



class PhotonJetAna : public edm::EDAnalyzer {
   public:
      explicit PhotonJetAna(const edm::ParameterSet&);
      ~PhotonJetAna();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

       bool    findGoodVtx; 
       math::XYZPoint vtx_;
       int          pdgId_;            // PDG ID of expected MC particle
       std::vector<int> otherPdgIds_;
 
  protected:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void storeGeneral(const edm::Event&);
      virtual void storeVertex(const edm::Event&); 
      virtual void storeJets(const edm::Event&); 
      virtual void YJstoreJets(const edm::Event&); 
      
      virtual int selectStorePhotons(const edm::Event&,const edm::EventSetup&);
      virtual int storePhotons(const edm::Event&,const edm::EventSetup&  ,pat::PhotonCollection &);

      virtual Float_t getGenCalIso(edm::Handle<reco::GenParticleCollection> handle,
                              reco::GenParticleCollection::const_iterator thisPho, const Float_t dRMax=0.4, bool removeMu=true, bool removeNu=false);

 virtual Float_t getGenTrkIso(edm::Handle<reco::GenParticleCollection> handle,
                              reco::GenParticleCollection::const_iterator thisPho, const Float_t dRMax=0.4);

      // ----------member data ---------------------------
      edm::InputTag VertexProducerC;      // vertecies producer
      edm::InputTag BeamSpotProducerC;  
      edm::InputTag JetProducerC;  
      edm::InputTag triggerEventProducerC;  
      edm::InputTag genParticlesProducerC;  
      edm::InputTag YJJetProducerC; 
      edm::InputTag PhotonProducerC; 
      edm::InputTag ebReducedRecHitCollectionC;
      edm::InputTag eeReducedRecHitCollectionC; 
      double ptMinC; 
      double etaMaxC;
      double hadEmMaxC;
      bool         muInGenCaloIso_; 

 bool isMC_;


      Bool_t doGenParticles_;

 
      
      TTree    *tree_;
      
      //general 
      Int_t    run_;
      Int_t    event_;
      Int_t    orbit_;
      Int_t    bx_;
      Int_t    lumis_;
      Bool_t   isData_;

      //beam spot    
      Float_t  bspotPos_[3];

      //vertex

      Int_t    nVtxGood_;   
      Int_t    nVtxNotFake_;
      Float_t  vertexX_[50];
      Float_t  vertexY_[50];
      Float_t  vertexZ_[50]; 
      Float_t  vertexXError_[50];
      Float_t  vertexYError_[50];
      Float_t  vertexZError_[50]; 
      Float_t  vertexChi2_[50];
      Float_t  vertexNormChi2_[50];
      Float_t  vertexNdof_[50]; 
  
      Int_t    vertexNTrk_[50];
      Int_t    vertexNTrkWeight05_[50];


      //Jet

      Int_t    nJet_;
      Int_t    jetTrg_[maxP][14];
      Int_t    jetAlgo_[maxP];
      Float_t  jetEn_[maxP];
      Float_t  jetPt_[maxP];
      Float_t  jetEta_[maxP];
      Float_t  jetPhi_[maxP];
      Float_t  jetEt_[maxP];
      Float_t  jetRawPt_[maxP];
      Float_t  jetRawEn_[maxP];
      Float_t  jetCharge_[maxP];
      Float_t  jetCHF_[maxP];
      Float_t  jetNHF_[maxP];
      Float_t  jetCEF_[maxP];
      Float_t  jetNEF_[maxP];
      Int_t    jetNCH_[maxP];
      Float_t  jetHFHAE_[maxP];
      Float_t  jetHFEME_[maxP];
      Int_t    jetPartonID_[maxP];
      Int_t    jetNConstituents_[maxP];
      Float_t  jetTrackCountHiEffBJetTags_[maxP];
      Float_t  jetTrackCountHiPurBJetTags_[maxP];
      Float_t  jetSimpleSVHiEffBJetTags_[maxP];
      Float_t  jetSimpleSVHiPurBJetTags_[maxP];
      
      Int_t    jetGenJetIndex_[maxP];
      Float_t  jetGenJetEn_[maxP];
      Float_t  jetGenJetPt_[maxP];
      Float_t  jetGenJetEta_[maxP];
      Float_t  jetGenJetPhi_[maxP];
      Int_t    jetGenPartonID_[maxP];
      Float_t  jetGenEn_[maxP];
      Float_t  jetGenPt_[maxP];
      Float_t  jetGenEta_[maxP];
      Float_t  jetGenPhi_[maxP];

      //YJJet

      Int_t    nPatJet_;
      Float_t  PatJetEn_[maxP];
      Float_t  PatJetPt_[maxP];
      Float_t  PatJetEta_[maxP];
      Float_t  PatJetPhi_[maxP];
      Float_t  PatJetEt_[maxP];
      Float_t  PatJetRawPt_[maxP];
      Float_t  PatJetRawEn_[maxP];
      Float_t  PatJetCharge_[maxP];
      Float_t  PatJetCHF_[maxP];
      Float_t  PatJetNHF_[maxP];
      Float_t  PatJetCEF_[maxP];
      Float_t  PatJetNEF_[maxP];
      Int_t    PatJetNCH_[maxP];
      Int_t    PatJetPartonID_[maxP];

      Int_t    PatJetGenJetIndex_[maxP];
      Float_t  PatJetGenJetEn_[maxP];
      Float_t  PatJetGenJetPt_[maxP];
      Float_t  PatJetGenJetEta_[maxP];
      Float_t  PatJetGenJetPhi_[maxP];
      Int_t    PatJetGenPartonID_[maxP];
      Float_t  PatJetGenEn_[maxP];
      Float_t  PatJetGenPt_[maxP];
      Float_t  PatJetGenEta_[maxP];
      Float_t  PatJetGenPhi_[maxP];



      Int_t    nPFJet_;
      Float_t  PFJetEn_[maxP];
      Float_t  PFJetPt_[maxP];
      Float_t  PFJetEta_[maxP];
      Float_t  PFJetPhi_[maxP];
      Float_t  PFJetEt_[maxP];
      Float_t  PFJetRawPt_[maxP];
      Float_t  PFJetRawEn_[maxP];
      Float_t  PFJetCharge_[maxP];
      Float_t  PFJetCHF_[maxP];
      Float_t  PFJetNHF_[maxP];
      Float_t  PFJetCEF_[maxP];
      Float_t  PFJetNEF_[maxP];
      Int_t    PFJetNCH_[maxP];
   
      //photon

      Int_t    nPho_;
      Float_t  PhoP_[maxP]; 
      Float_t  PhoEt_[maxP]; 
      Float_t  PhoEnergy_[maxP]; 
      Float_t  PhoPx_[maxP]; 
      Float_t  PhoPy_[maxP]; 
      Float_t  PhoPz_[maxP]; 
      Float_t  PhoPt_[maxP]; 
      Float_t  PhoEta_[maxP]; 
      Float_t  PhoPhi_[maxP]; 
      Float_t  PhoR9_[maxP]; 
      Float_t  PhoPhiWidth_[maxP]; 
      Float_t  PhoEtaWidth_[maxP]; 
      Float_t  PhoScPhi_[maxP]; 
      Float_t  PhoScEta_[maxP]; 
      Float_t  PhoESRatio_[maxP];
      Float_t  PhoSigmaIetaIeta_[maxP];
      Float_t  PhoSigmaIphiIphi_[maxP];
      Float_t  PhoSeedTime_[maxP];
      Float_t  PhoE2overe9_[maxP];
      Float_t  PhohadronicOverEm_[maxP];
      Float_t  PhoecalRecHitSumEtConeDR04_[maxP];
      Float_t  PhohcalTowerSumEtConeDR04_[maxP];
      Float_t  PhotrkSumPtHollowConeDR04_[maxP]; 
      Float_t  PhoisConverted_[maxP]; 
      Float_t  PhohasPixelSeed_[maxP]; 
      Float_t  PhoisGenMatched_[maxP];
      
      Float_t  PhogenMomId_[maxP]; 
      Float_t  PhogenGrandMomId_[maxP];
      Float_t  PhogenNSiblings_[maxP];
      Float_t  PhogenMatchedPt_[maxP];
      Float_t  PhogenMatchedEta_[maxP];
      Float_t  PhogenMatchedPhi_[maxP];
      Float_t  PhogenCalIsoDR04_[maxP];
      Float_t  PhogenTrkIsoDR04_[maxP];
      Float_t  PhogenIsoDR04_[maxP];
 
      /*
      Float_t  pdf_[7];
      Float_t  pthat_;
      Float_t  processID_;
      Int_t    nHLT_;
      Int_t    HLT_[maxP];
      Int_t    HLTIndex_[50];
      Float_t  bspotPos_[3];
      Int_t    nVtx_;
      Float_t  vtx_[50][3];
      Int_t    vtxNTrk_[50];
      Int_t    vtxNDF_[50];
      Float_t  vtxD0_[50];
      Int_t    IsVtxGood_;
      Int_t    nVtxBS_;
      Float_t  vtxbs_[50][3];
      Float_t  vtxbsPtMod_[50];
      Float_t  vtxbsSumPt2_[50];
*/

};

