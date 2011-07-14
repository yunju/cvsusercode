// -*- C++ -*-
//
// Package:    PhotonJetAna
// Class:      PhotonJetAna
// 
/**\class PhotonJetAna PhotonJetAna.cc YJAnaPhotonJet/PhotonJetAna/src/PhotonJetAna.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yun-Ju Lu,27 2-004,+41227676186,
//         Created:  Fri Jul  8 15:56:55 CEST 2011
// $Id: PhotonJetAna.cc,v 1.1 2011/07/12 17:06:56 yunju Exp $
//
//


// system include files
#include "YJAnaPhotonJet/PhotonJetAna/plugins/PhotonJetAna.h"
 
 #include "DataFormats/HepMCCandidate/interface/GenParticle.h"
 #include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
 #include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"

  #include "FWCore/Framework/interface/EDAnalyzer.h"

 //For pixel RecHit handle                                                                                                      
  #include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
  #include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
  #include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
  #include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

  #include "CommonTools/Utils/interface/PtComparator.h"


//ROOT includes
 #include <Math/VectorUtil.h>
 #include <TLorentzVector.h>

using namespace ROOT::Math::VectorUtil;


PhotonJetAna::PhotonJetAna(const edm::ParameterSet& ps)

{
   //now do what ever initialization is needed
    //c inputtag to cfi   
     VertexProducerC       = ps.getParameter<edm::InputTag>("VertexProducerPY");
     BeamSpotProducerC     = ps.getParameter<edm::InputTag>("BeamSpotProducerPY"); 
     JetProducerC          = ps.getParameter<edm::InputTag>("JetProducerPY");   
     triggerEventProducerC = ps.getParameter<edm::InputTag>("triggerEventProducerPY");   
     genParticlesProducerC = ps.getParameter<edm::InputTag>("genParticlesProducerPY");   
     YJJetProducerC        = ps.getParameter<edm::InputTag>("YJJetProducerPY") ;
     PhotonProducerC       = ps.getParameter<edm::InputTag>("PhotonProducerPY") ;   
     ebReducedRecHitCollectionC = ps.getParameter<edm::InputTag>("ebReducedRecHitCollectionPY") ;
     eeReducedRecHitCollectionC = ps.getParameter<edm::InputTag>("eeReducedRecHitCollectionPY") ;
     ptMinC                = ps.getUntrackedParameter<double>("GammaPtMin", 15);
     etaMaxC               = ps.getUntrackedParameter<double>("GammaEtaMax",3);
     hadEmMaxC             = ps.getUntrackedParameter<double>("GammaHadEmMax",0.05);
     pdgId_                           = ps.getUntrackedParameter<int>("pdgId", 22);

     doGenParticles_   = ps.getParameter<bool>("doGenParticles");
 muInGenCaloIso_                  = ps.getUntrackedParameter<bool>("muInGenCaloIso", 10);   
   
 isMC_=true;

   Service<TFileService> fs;
    tree_ = fs->make<TTree>("EventTree", "Event data");
    tree_->Branch("run", &run_, "run/I");
   tree_->Branch("event", &event_, "event/I");
   tree_->Branch("orbit", &orbit_, "orbit/I");
   tree_->Branch("bx", &bx_, "bx/I");
   tree_->Branch("lumis", &lumis_, "lumis/I");
   tree_->Branch("isData", &isData_, "isData/O");

   //beam spot

   tree_->Branch("bspotPos", bspotPos_, "bspotPos[3]/F");

   //vertex

   tree_->Branch("nVtxGood", &nVtxGood_, "nVtxGood/I");
   tree_->Branch("nVtxNotFake", &nVtxNotFake_, "nVtxNotFake/I");
   tree_->Branch("vertexX", vertexX_, "vertexX[nVtxNotFake]/F");
   tree_->Branch("vertexY", vertexY_, "vertexY[nVtxNotFake]/F");
   tree_->Branch("vertexZ", vertexZ_, "vertexZ[nVtxNotFake]/F");
   tree_->Branch("vertexXError", vertexXError_, "vertexXError[nVtxNotFake]/F");
   tree_->Branch("vertexYError", vertexYError_, "vertexYError[nVtxNotFake]/F");
   tree_->Branch("vertexZError", vertexZError_, "vertexZError[nVtxNotFake]/F");
   tree_->Branch("vertexChi2", vertexChi2_, "vertexChi2[nVtxNotFake]/F");
   tree_->Branch("vertexNormChi2", vertexNormChi2_, "vertexNormChi2[nVtxNotFake]/F");
   tree_->Branch("vertexNdof", vertexNdof_, "vertexNdof[nVtxNotFake]/F");
   tree_->Branch("vertexNTrk", vertexNTrk_, "vertexNTrk[nVtxNotFake]/F");
   tree_->Branch("vertexNTrkWeight05", vertexNTrkWeight05_, "vertexNTrkWeight05[nVtxNotFake]/F");

   //Jets
/*
   tree_->Branch("nJet", &nJet_, "nJet/I");
   tree_->Branch("jetTrg", jetTrg_, "jetTrg[nJet][14]/I");
   tree_->Branch("jetEn", jetEn_, "jetEn[nJet]/F");
   tree_->Branch("jetPt", jetPt_, "jetPt[nJet]/F");
   tree_->Branch("jetEta", jetEta_, "jetEta[nJet]/F");
   tree_->Branch("jetPhi", jetPhi_, "jetPhi[nJet]/F");
   tree_->Branch("jetEt", jetEt_, "jetEt[nJet]/F");
   tree_->Branch("jetRawPt", jetRawPt_, "jetRawPt[nJet]/F");
   tree_->Branch("jetRawEn", jetRawEn_, "jetRawEn[nJet]/F");
   tree_->Branch("jetCHF", jetCHF_, "jetCHF[nJet]/F");
   tree_->Branch("jetNHF", jetNHF_, "jetNHF[nJet]/F");
   tree_->Branch("jetCEF", jetCEF_, "jetCEF[nJet]/F");
   tree_->Branch("jetNEF", jetNEF_, "jetNEF[nJet]/F");
   tree_->Branch("jetNCH", jetNCH_, "jetNCH[nJet]/I");
   tree_->Branch("jetHFHAE", jetHFHAE_, "jetHFHAE[nJet]/F");
   tree_->Branch("jetHFEME", jetHFEME_, "jetHFEME[nJet]/F");
   tree_->Branch("jetNConstituents", jetNConstituents_, "jetNConstituents[nJet]/I");
   tree_->Branch("jetTrackCountHiEffBJetTags", jetTrackCountHiEffBJetTags_, "jetTrackCountHiEffBJetTags[nJet]/F");
   tree_->Branch("jetTrackCountHiPurBJetTags", jetTrackCountHiPurBJetTags_, "jetTrackCountHiPurBJetTags[nJet]/F");
   tree_->Branch("jetSimpleSVHiEffBJetTags", jetSimpleSVHiEffBJetTags_, "jetSimpleSVHiEffBJetTags[nJet]/F");
   tree_->Branch("jetSimpleSVHiPurBJetTags", jetSimpleSVHiPurBJetTags_, "jetSimpleSVHiPurBJetTags[nJet]/F");
   if (doGenParticles_) {
     tree_->Branch("jetPartonID", jetPartonID_, "jetPartonID[nJet]/I");
     tree_->Branch("jetGenJetIndex", jetGenJetIndex_, "jetGenJetIndex[nJet]/I");
     tree_->Branch("jetGenJetEn", jetGenJetEn_, "jetGenJetEn[nJet]/F");
     tree_->Branch("jetGenJetPt", jetGenJetPt_, "jetGenJetPt[nJet]/F");
     tree_->Branch("jetGenJetEta", jetGenJetEta_, "jetGenJetEta[nJet]/F");
     tree_->Branch("jetGenJetPhi", jetGenJetPhi_, "jetGenJetPhi[nJet]/F");
     tree_->Branch("jetGenPartonID", jetGenPartonID_, "jetGenPartonID[nJet]/I");
     tree_->Branch("jetGenEn", jetGenEn_, "jetGenEn[nJet]/F");
     tree_->Branch("jetGenPt", jetGenPt_, "jetGenPt[nJet]/F");
     tree_->Branch("jetGenEta", jetGenEta_, "jetGenEta[nJet]/F");
     tree_->Branch("jetGenPhi", jetGenPhi_, "jetGenPhi[nJet]/F");
   }
*/
//YJJet

     tree_->Branch("nPatJet", &nPatJet_, "nPatJet/I");
     tree_->Branch("PatJetEn", PatJetEn_, "PatJetEn[nPatJet]/F");
     tree_->Branch("PatJetPt", PatJetPt_, "PatJetPt[nPatJet]/F");
     tree_->Branch("PatJetEta", PatJetEta_, "PatJetEta[nPatJet]/F");
     tree_->Branch("PatJetPhi", PatJetPhi_, "PatJetPhi[nPatJet]/F");
     tree_->Branch("PatJetEt", PatJetEt_, "PatJetEt[nPatJet]/F");
     tree_->Branch("PatJetRawPt", PatJetRawPt_, "PatJetRawPt[nPatJet]/F");
     tree_->Branch("PatJetRawEn", PatJetRawEn_, "PatJetRawEn[nPatJet]/F");
     tree_->Branch("PatJetCHF", PatJetCHF_, "PatJetCHF[nPatJet]/F");
     tree_->Branch("PatJetNHF", PatJetNHF_, "PatJetNHF[nPatJet]/F");
     tree_->Branch("PatJetCEF", PatJetCEF_, "PatJetCEF[nPatJet]/F");
     tree_->Branch("PatJetNEF", PatJetNEF_, "PatJetNEF[nPatJet]/F");
     tree_->Branch("PatJetNCH", PatJetNCH_, "PatJetNCH[nPatJet]/I");
     tree_->Branch("PatJetPartonID", PatJetPartonID_, "PatJetPartonID[nPatJet]/I");

     tree_->Branch("PatJetGenJetIndex", PatJetGenJetIndex_, "PatJetGenJetIndex[nPatJet]/I");
     tree_->Branch("PatJetGenJetEn", PatJetGenJetEn_, "PatJetGenJetEn[nPatJet]/F");
     tree_->Branch("PatJetGenJetPt", PatJetGenJetPt_, "PatJetGenJetPt[nPatJet]/F");
     tree_->Branch("PatJetGenJetEta", PatJetGenJetEta_, "PatJetGenJetEta[nPatJet]/F");
     tree_->Branch("PatJetGenJetPhi", PatJetGenJetPhi_, "PatJetGenJetPhi[nPatJet]/F");
     tree_->Branch("PatJetGenPartonID", PatJetGenPartonID_, "PatJetGenPartonID[nPatJet]/I");
     tree_->Branch("PatJetGenEn", PatJetGenEn_, "PatJetGenEn[nPatJet]/F");
     tree_->Branch("PatJetGenPt", PatJetGenPt_, "PatJetGenPt[nPatJet]/F");
     tree_->Branch("PatJetGenEta", PatJetGenEta_, "PatJetGenEta[nPatJet]/F");
     tree_->Branch("PatJetGenPhi", PatJetGenPhi_, "PatJetGenPhi[nPatJet]/F");
/*
     tree_->Branch("nPFJet", &nPFJet_, "nPFJet/I");
     tree_->Branch("PFJetEn", PFJetEn_, "PFJetEn[nPFJet]/F");
     tree_->Branch("PFJetPt", PFJetPt_, "PFJetPt[nPFJet]/F");
     tree_->Branch("PFJetEta", PFJetEta_, "PFJetEta[nPFJet]/F");
     tree_->Branch("PFJetPhi", PFJetPhi_, "PFJetPhi[nPFJet]/F");
     tree_->Branch("PFJetEt", PFJetEt_, "PFJetEt[nPFJet]/F");
     tree_->Branch("PFJetRawPt", PFJetRawPt_, "PFJetRawPt[nPFJet]/F");
     tree_->Branch("PFJetRawEn", PFJetRawEn_, "PFJetRawEn[nPFJet]/F");
     tree_->Branch("PFJetCHF", PFJetCHF_, "PFJetCHF[nPFJet]/F");
     tree_->Branch("PFJetNHF", PFJetNHF_, "PFJetNHF[nPFJet]/F");
     tree_->Branch("PFJetCEF", PFJetCEF_, "PFJetCEF[nPFJet]/F");
     tree_->Branch("PFJetNEF", PFJetNEF_, "PFJetNEF[nPFJet]/F");
     tree_->Branch("PFJetNCH", PFJetNCH_, "PFJetNCH[nPFJet]/I");
*/

     //photon

     tree_->Branch("nPho", &nPho_, "nPho/I");
     tree_->Branch("PhoP", PhoP_, "PhoP[nPho]/F");
     tree_->Branch("PhoEt", PhoEt_, "PhoEt[nPho]/F");
     tree_->Branch("PhoEnergy", PhoEnergy_, "PhoEnergy[nPho]/F");
     tree_->Branch("PhoPx", PhoPx_, "PhoPx[nPho]/F");
     tree_->Branch("PhoPy", PhoPy_, "PhoPy[nPho]/F");
     tree_->Branch("PhoPz", PhoPz_, "PhoPz[nPho]/F");
     tree_->Branch("PhoEta", PhoEta_, "PhoEta[nPho]/F");
     tree_->Branch("PhoPhi", PhoPhi_, "PhoPhi[nPho]/F");
     tree_->Branch("PhoR9", PhoR9_, "PhoR9[nPho]/F");
     tree_->Branch("PhoPhiWidth", PhoPhiWidth_, "PhoPhiWidth[nPho]/F");
     tree_->Branch("PhoEtaWidth", PhoEtaWidth_, "PhoEtaWidth[nPho]/F");
     tree_->Branch("PhoScPhi", PhoScPhi_, "PhoScPhi[nPho]/F");
     tree_->Branch("PhoScEta", PhoScEta_, "PhoScEta[nPho]/F");
     tree_->Branch("PhoESRatio", PhoESRatio_, "PhoESRatio[nPho]/F");
    
     tree_->Branch("PhoSigmaIetaIeta", PhoSigmaIetaIeta_, "PhoSigmaIetaIeta[nPho]/F");
     tree_->Branch("PhoSigmaIphiIphi", PhoSigmaIphiIphi_, "PhoSigmaIphiIphi[nPho]/F");
     tree_->Branch("PhoSeedTime", PhoSeedTime_, "PhoSeedTime[nPho]/F");
     tree_->Branch("PhoE2overe9", PhoE2overe9_, "PhoE2overe9[nPho]/F");
     tree_->Branch("PhohadronicOverEm", PhohadronicOverEm_, "PhohadronicOverEm[nPho]/F");
     tree_->Branch("PhoecalRecHitSumEtConeDR04", PhoecalRecHitSumEtConeDR04_, "PhoecalRecHitSumEtConeDR04[nPho]/F");
     tree_->Branch("PhohcalTowerSumEtConeDR04", PhohcalTowerSumEtConeDR04_, "PhohcalTowerSumEtConeDR04[nPho]/F");
     tree_->Branch("PhotrkSumPtHollowConeDR04", PhotrkSumPtHollowConeDR04_, "PhotrkSumPtHollowConeDR04[nPho]/F");
     tree_->Branch("PhoisConverted", PhoisConverted_, "PhoisConverted[nPho]/F");
     tree_->Branch("PhohasPixelSeed", PhohasPixelSeed_, "PhohasPixelSeed[nPho]/F");
    
     tree_->Branch("PhoisGenMatched", PhoisGenMatched_, "PhoisGenMatched[nPho]/F");

     tree_->Branch("PhogenMomId", PhogenMomId_, "PhogenMomId[nPho]/F");
     tree_->Branch("PhogenGrandMomId", PhogenGrandMomId_, "PhogenGrandMomId[nPho]/F");
     tree_->Branch("PhogenNSiblings", PhogenNSiblings_, "PhogenNSiblings[nPho]/F");
     tree_->Branch("PhogenMatchedPt", PhogenMatchedPt_, "PhogenMatchedPt[nPho]/F");
     tree_->Branch("PhogenMatchedEta", PhogenMatchedEta_, "PhogenMatchedEta[nPho]/F");
     tree_->Branch("PhogenMatchedPhi", PhogenMatchedPhi_, "PhogenMatchedPhi[nPho]/F");
     tree_->Branch("PhogenCalIsoDR04", PhogenCalIsoDR04_, "PhogenCalIsoDR04[nPho]/F");
     tree_->Branch("PhogenTrkIsoDR04", PhogenTrkIsoDR04_, "PhogenTrkIsoDR04[nPho]/F");
     tree_->Branch("PhogenIsoDR04", PhogenIsoDR04_, "PhogenIsoDR04[nPho]/F");


     





 
  
    

     




























/*

   tree_->Branch("nHLT", &nHLT_, "nHLT/I");
   tree_->Branch("HLT", HLT_, "HLT[nHLT]/I");
   tree_->Branch("HLTIndex", HLTIndex_, "HLTIndex[50]/I");
   tree_->Branch("bspotPos", bspotPos_, "bspotPos[3]/F");
   tree_->Branch("nVtx", &nVtx_, "nVtx/I");
   tree_->Branch("vtx", vtx_, "vtx[nVtx][3]/F");
   tree_->Branch("vtxNTrk", vtxNTrk_, "vtxNTrk[nVtx]/I");
   tree_->Branch("vtxNDF", vtxNDF_, "vtxNDF[nVtx]/I");
   tree_->Branch("vtxD0", vtxD0_, "vtxD0[nVtx]/F");
   tree_->Branch("IsVtxGood", &IsVtxGood_, "IsVtxGood/I");
   tree_->Branch("nVtxBS", &nVtxBS_, "nVtxBS/I");
   tree_->Branch("vtxbs", vtxbs_, "vtxbs[nVtxBS][3]/F");
   tree_->Branch("vtxbsPtMod", vtxbsPtMod_, "vtxbsPtMod[nVtxBS]/F");
*/


}


PhotonJetAna::~PhotonJetAna()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhotonJetAna::analyze(const edm::Event& e, const edm::EventSetup& iSetup)
{

     
   using namespace edm;
   storeGeneral(e);
   storeVertex(e);
  // storeJets(e);
   YJstoreJets(e);
   bool foundPhotons = selectStorePhotons(e,iSetup);
  // if (foundPhotons){
    tree_->Fill(); 
  // }
}

void PhotonJetAna::storeGeneral(const edm::Event& e)
{
   run_    = e.id().run();
   event_  = e.id().event();
   orbit_  = e.orbitNumber();
   bx_     = e.bunchCrossing();
   lumis_  = e.luminosityBlock();
   isData_ = e.isRealData();


}

void PhotonJetAna::YJstoreJets(const edm::Event& e)
{
   Handle<View<pat::Jet> > colojets;
   e.getByLabel(YJJetProducerC,colojets);

 
   nPatJet_ = 0;
   for (View<pat::Jet>::const_iterator iJet = colojets->begin(); iJet != colojets->end(); ++iJet) 
   {
     //cout<<"OK?"<<endl;

      PatJetEn_[nPatJet_]     = iJet->energy();
      PatJetPt_[nPatJet_]     = iJet->pt();
      PatJetEta_[nPatJet_]    = iJet->eta();
      PatJetPhi_[nPatJet_]    = iJet->phi();
      PatJetCharge_[nPatJet_] = iJet->jetCharge();
      PatJetEt_[nPatJet_]     = iJet->et();
     //  cout<<"OK2?"<<endl;
/*
      PatJetRawPt_[nPatJet_]  = (*iJet).correctedJet("Uncorrected").pt();
      PatJetRawEn_[nPatJet_]  = (*iJet).correctedJet("Uncorrected").energy();
*/


      PatJetCEF_[nPatJet_]    = iJet->chargedEmEnergyFraction();

      PatJetNEF_[nPatJet_]    = iJet->neutralEmEnergyFraction();
      PatJetCHF_[nPatJet_]    = iJet->chargedHadronEnergyFraction();
      PatJetNHF_[nPatJet_]    = iJet->neutralHadronEnergyFraction();
      PatJetNCH_[nPatJet_]    = iJet->chargedMultiplicity();
      PatJetPartonID_[nPatJet_] = iJet->partonFlavour();
 
      edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle; 
      e.getByLabel(genParticlesProducerC,genParticlesHandle);

    

 
      if ( genParticlesHandle.isValid()&&(*iJet).genParton()) 
      {
          PatJetGenPartonID_[nPatJet_] = (*iJet).genParton()->pdgId();
          PatJetGenEn_[nPatJet_]       = (*iJet).genParton()->energy();
          PatJetGenPt_[nPatJet_]       = (*iJet).genParton()->pt();
          PatJetGenEta_[nPatJet_]      = (*iJet).genParton()->eta();
          PatJetGenPhi_[nPatJet_]      = (*iJet).genParton()->phi();


      }
      

      if (genParticlesHandle.isValid()&&(*iJet).genJet()) 
      {
          PatJetGenJetIndex_[nPatJet_] = 1;
	  PatJetGenJetEn_[nPatJet_] = (*iJet).genJet()->energy();
	  PatJetGenJetPt_[nPatJet_] = (*iJet).genJet()->pt();
	  PatJetGenJetEta_[nPatJet_] = (*iJet).genJet()->eta();
	  PatJetGenJetPhi_[nPatJet_] = (*iJet).genJet()->phi();

      }
       


     nPatJet_++;
   }


/*      
   Handle<JetCollection> pfjets;
   e.getByLabel("ak5PFJets",pfjets);
   
   Handle<JetCollection> genjets;
   e.getByLabel("ak5GenJets",genjets);
 */  


}


void PhotonJetAna::storeJets(const edm::Event& e)
{
   edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle; 
   e.getByLabel(genParticlesProducerC,genParticlesHandle);
/*
   edm::Handle<TriggerEvent> triggerEventHandle; 
   e.getByLabel(triggerEventProducerC,triggerEventHandle);
 
  const TriggerObjectMatch *jetTriggerMatch1(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet15U"));
  const TriggerObjectMatch *jetTriggerMatch2(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet30U"));
  const TriggerObjectMatch *jetTriggerMatch3(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet50U"));
  const TriggerObjectMatch *jetTriggerMatch4(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet70U"));
  const TriggerObjectMatch *jetTriggerMatch5(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet70Uv2"));
  const TriggerObjectMatch *jetTriggerMatch6(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet100U"));
  const TriggerObjectMatch *jetTriggerMatch7(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet100Uv2"));
  const TriggerObjectMatch *jetTriggerMatch8(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet140Uv1"));
  const TriggerObjectMatch *jetTriggerMatch9(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet15Uv3"));
  const TriggerObjectMatch *jetTriggerMatch10(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet30Uv3"));
  const TriggerObjectMatch *jetTriggerMatch11(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet50Uv3"));
  const TriggerObjectMatch *jetTriggerMatch12(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet70Uv3"));
  const TriggerObjectMatch *jetTriggerMatch13(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet100Uv3"));
  const TriggerObjectMatch *jetTriggerMatch14(triggerEventHandle->triggerObjectMatchResult("jetTriggerMatchHLTJet140Uv3"));

*/
   nJet_ = 0;
   const TriggerMatchHelper matchHelper;
   
    edm::Handle<View<pat::Jet> > jetHandle;
    e.getByLabel(JetProducerC,jetHandle);

  if (jetHandle.isValid())
   

     cout<<"Where"<<endl; 
    for (View<pat::Jet>::const_iterator iJet = jetHandle->begin(); iJet != jetHandle->end(); ++iJet) {
      
      if (iJet->pt() < 15) continue;

      edm::RefToBase<pat::Jet> jetRef = jetHandle->refAt(nJet_);
      reco::CandidateBaseRef jetBaseRef(jetRef);
     /* 
      const TriggerObjectRef jetTrigRef1( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch1, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef2( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch2, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef3( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch3, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef4( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch4, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef5( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch5, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef6( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch6, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef7( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch7, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef8( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch8, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef9( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch9, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef10( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch10, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef11( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch11, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef12( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch12, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef13( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch13, e, *triggerEventHandle ) );
      const TriggerObjectRef jetTrigRef14( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch14, e, *triggerEventHandle ) );
      jetTrg_[nJet_][0] = (jetTrigRef1.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][1] = (jetTrigRef2.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][2] = (jetTrigRef3.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][3] = (jetTrigRef4.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][4] = (jetTrigRef5.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][5] = (jetTrigRef6.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][6] = (jetTrigRef7.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][7] = (jetTrigRef8.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][8] = (jetTrigRef9.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][9] = (jetTrigRef10.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][10] = (jetTrigRef11.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][11] = (jetTrigRef12.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][12] = (jetTrigRef13.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][13] = (jetTrigRef14.isAvailable()) ? 1 : -99;
*/
      jetEn_[nJet_]     = iJet->energy();
      jetPt_[nJet_]     = iJet->pt();
      jetEta_[nJet_]    = iJet->eta();
      jetPhi_[nJet_]    = iJet->phi();
      jetCharge_[nJet_] = iJet->jetCharge();
      jetEt_[nJet_]     = iJet->et();
      jetRawPt_[nJet_]  = (*iJet).correctedJet("Uncorrected").pt();
      jetRawEn_[nJet_]  = (*iJet).correctedJet("Uncorrected").energy();

      jetCEF_[nJet_]    = iJet->chargedEmEnergyFraction();
      jetNEF_[nJet_]    = iJet->neutralEmEnergyFraction();
      jetCHF_[nJet_]    = iJet->chargedHadronEnergyFraction();
      jetNHF_[nJet_]    = iJet->neutralHadronEnergyFraction();
      jetHFHAE_[nJet_]  = iJet->HFHadronEnergy();
      jetHFEME_[nJet_]  = iJet->HFEMEnergy();
      jetNCH_[nJet_]    = iJet->chargedMultiplicity();
      jetNConstituents_[nJet_] = iJet->getPFConstituents().size();

      // b-tagging
      jetTrackCountHiEffBJetTags_[nJet_] = iJet->bDiscriminator("trackCountingHighEffBJetTags");
      jetTrackCountHiPurBJetTags_[nJet_] = iJet->bDiscriminator("trackCountingHighPurBJetTags");
      jetSimpleSVHiEffBJetTags_[nJet_]   = iJet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      jetSimpleSVHiPurBJetTags_[nJet_]   = iJet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
      
      // gen jet and parton
      jetPartonID_[nJet_] = iJet->partonFlavour();

      jetGenPartonID_[nJet_]    = -99;
      if (!isData_ && genParticlesHandle.isValid() ) {
        if ((*iJet).genParton()) {
          jetGenPartonID_[nJet_] = (*iJet).genParton()->pdgId();
          jetGenEn_[nJet_]       = (*iJet).genParton()->energy();
          jetGenPt_[nJet_]       = (*iJet).genParton()->pt();
          jetGenEta_[nJet_]      = (*iJet).genParton()->eta();
          jetGenPhi_[nJet_]      = (*iJet).genParton()->phi();
	}
      }

      jetGenJetIndex_[nJet_] = -1;
      jetGenJetEn_[nJet_]    = -1;
      jetGenJetPt_[nJet_]    = -999;
      jetGenJetEta_[nJet_]   = -999;
      jetGenJetPhi_[nJet_]   = -999;
      
      if (!isData_ && genParticlesHandle.isValid() ) {
	if ((*iJet).genJet()) {
	  jetGenJetIndex_[nJet_] = 1;
	  jetGenJetEn_[nJet_] = (*iJet).genJet()->energy();
	  jetGenJetPt_[nJet_] = (*iJet).genJet()->pt();
	  jetGenJetEta_[nJet_] = (*iJet).genJet()->eta();
	  jetGenJetPhi_[nJet_] = (*iJet).genJet()->phi();
	}
      }
      nJet_++;
    }

  
}




void PhotonJetAna::storeVertex(const edm::Event& e){
	///////////////////////////////////////////////////////////////////////
	// Vertex Section: store BeamSpot and Primary Vertex of the event    //
	///////////////////////////////////////////////////////////////////////
  nVtxNotFake_ = 0;
	
//  cout <<" start storeVertex "<<endl;
  findGoodVtx          = false;
	// Get the Beam Spot
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  e.getByLabel(BeamSpotProducerC,recoBeamSpotHandle);
  beamSpot = *recoBeamSpotHandle;
 
  // beam spot position
  bspotPos_[0] = beamSpot.x0();
  bspotPos_[1] = beamSpot.y0();
  bspotPos_[2] = beamSpot.z0();

 
	// Get the primary event vertex
  Handle<reco::VertexCollection> vertexHandle;
  e.getByLabel(VertexProducerC, vertexHandle);
  reco::VertexCollection vertexCollection = *(vertexHandle.product());
  vtx_.SetXYZ(0.,0.,0.);
  double chi2(-1), ndof(-1), normChi2(-1), vtxXError(-1),  vtxYError(-1), vtxZError(-1);
  int vtxNTrk(0), vtxNTrkWeight05(0), nVtxGood(0);


   

  if (vertexCollection.size()>0) 
  { 
  
    reco::VertexCollection::const_iterator vert = vertexCollection.begin(); 

    for( ; vert!=vertexCollection.end(); ++vert)
    {
       Double_t vertex_rho = sqrt(vert->x()*vert->x()+
				 vert->y()*vert->y());
 
      if ( !(vert->isFake()) && vert->ndof() > 4 && fabs(vert->z()) < 24.0   && vertex_rho < 2.0 )
      {
	  if(nVtxGood==0)vtx_ = vert->position();  // set vtx_ to be the first and best vertex
	  nVtxGood++;	
	  findGoodVtx=true;
      }//(!(vert->isFake()) && vert->ndof() > 4 && fabs(vert->z()) < 24.0   && vertex_rho < 2.0)
      nVtxGood_=nVtxGood;
      
      if ( !(vert->isFake()) ) 
      {
	vtxXError = vert->xError();
	vtxYError = vert->yError();
	vtxZError = vert->zError();
	chi2      = vert->chi2();  
	ndof      = vert->ndof();  
	normChi2  = vert->normalizedChi2();  
	vtxNTrk   = vert->tracksSize();
	
	vtxNTrkWeight05 = 0;
	reco::Vertex::trackRef_iterator ittrk;
	for(ittrk = vert->tracks_begin(); ittrk!= vert->tracks_end(); ++ittrk)
	  if ( vert->trackWeight(*ittrk) > 0.5 ) vtxNTrkWeight05++;
	
	vertexX_ [nVtxNotFake_] = vert->position().X();
	vertexY_ [nVtxNotFake_] = vert->position().Y();
	vertexZ_ [nVtxNotFake_] = vert->position().Z();
	vertexXError_ [nVtxNotFake_] = 	vtxXError;
	vertexYError_ [nVtxNotFake_] = 	vtxYError;
	vertexZError_ [nVtxNotFake_] = 	vtxZError;
	vertexChi2_ [nVtxNotFake_] = chi2;
	vertexNormChi2_ [nVtxNotFake_] = normChi2;
	vertexNdof_ [nVtxNotFake_] = ndof; 
	vertexNTrk_ [nVtxNotFake_] =  vtxNTrk;
	vertexNTrkWeight05_ [nVtxNotFake_] = vtxNTrkWeight05;

	nVtxNotFake_++;

      
     }// if ( !(vert->isFake()) )
    } // for( ; vert!=vertexCollection.end(); ++vert)
  }//if (vertexCollection.size()>0)


	
}


int PhotonJetAna::selectStorePhotons(const edm::Event& e,const edm::EventSetup& iSetup){

  /////////////////////////////////////////////////////////////////////////////
  // Photon Section: store kMaxPhotons in the events as an array in the tree //
  /////////////////////////////////////////////////////////////////////////////
  // Get photon details  
  Handle<PhotonCollection> photonshandle;
  e.getByLabel(PhotonProducerC, photonshandle);   
  
  // Sort photons according to pt
  PhotonCollection myphotons;
  for (PhotonCollection::const_iterator phoItr = photonshandle->begin(); phoItr != photonshandle->end(); ++phoItr) {  
    myphotons.push_back(*phoItr);
  }
  
  GreaterByPt<Photon> pTComparator_;
  std::sort(myphotons.begin(), myphotons.end(), pTComparator_);


  //cout<<"Here1 ?"<<endl;
  // Get the compl. Photons                                                                                                                                                         
  //  return storePhotons(e,iSetup,myphotons,prefx);
  return storePhotons(e,iSetup,myphotons);
}

// int MultiPhotonAnalyzer::storePhotons(const edm::Event& e,const edm::EventSetup& iSetup,PhotonCollection & myphotons, const char* prefx){
int PhotonJetAna::storePhotons(const edm::Event& e,const edm::EventSetup& iSetup,PhotonCollection &myphotons){
   
    

  // Tools to get cluster shapes

  edm::Handle<EcalRecHitCollection> EBReducedRecHits;
  e.getByLabel(ebReducedRecHitCollectionC, EBReducedRecHits);
  edm::Handle<EcalRecHitCollection> EEReducedRecHits;
  e.getByLabel(eeReducedRecHitCollectionC, EEReducedRecHits); 
  // get the channel status from the DB
  edm::ESHandle<EcalChannelStatus> chStatus;
  iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  
  // BasicCluster isolation varialbe calculators
/*
  CxCalculator CxC(e,iSetup, basicClusterBarrel_, basicClusterEndcap_);
  RxCalculator RxC(e,iSetup, hbhe_, hf_, ho_);
  TxyCalculator Txy(e,iSetup,trackProducer_);
  dRxyCalculator dRxy(e,iSetup,trackProducer_);
*/
  EcalClusterLazyTools lazyTool(e, iSetup, ebReducedRecHitCollectionC, eeReducedRecHitCollectionC );   
  
  // kmaxphotons

  const int kMaxPhotons = 30; 
  
   
    
    
  //// pixel hit handle ////
  int nPixel =0 ;
  const SiPixelRecHitCollection* rechits;
  Handle<SiPixelRecHitCollection> rchts;
  e.getByLabel("siPixelRecHits",rchts);
  if(rchts.isValid()){
    rechits = rchts.product();
    for (SiPixelRecHitCollection::const_iterator it = rechits->begin(); it!=rechits->end();it++)
      {
	SiPixelRecHitCollection::DetSet hits = *it;
        DetId detId = DetId(hits.detId());
	SiPixelRecHitCollection::const_iterator recHitMatch = rechits->find(detId);
	const SiPixelRecHitCollection::DetSet recHitRange = *recHitMatch;
        for ( SiPixelRecHitCollection::DetSet::const_iterator recHitIterator = recHitRange.begin();
              recHitIterator != recHitRange.end(); ++recHitIterator) {
	  // add selection if needed, now all hits.                                                                                                                                        
	  nPixel++;
        }
      }
  } // if the siPixelRecHit collection is saved in data
   
  nPho_=0;
  for (PhotonCollection::const_iterator phoItr = myphotons.begin(); phoItr != myphotons.end(); ++phoItr) {  
     //cout<<"nPho_"<< nPho_<<endl;
    
      if(phoItr->pt() < ptMinC || fabs(phoItr->p4().eta()) > etaMaxC || phoItr->hadronicOverEm() > hadEmMaxC) 
      continue;
    // Dump photon kinematics and AOD
    Photon photon = Photon(*phoItr);
    // NOTE: since CMSSW_3_1_x all photons are corrected to the primary vertex
    //       hence, Photon::setVertex() leaves photon object unchanged

    if(findGoodVtx)photon.setVertex(vtx_);
   // cout<<"Where_1"<<endl; 
    PhoP_[nPho_] =  photon.p();
    PhoEt_[nPho_] =  photon.et();
    PhoEnergy_[nPho_] =  photon.energy();
    PhoPx_[nPho_] =  photon.px();
    PhoPy_[nPho_] =  photon.py();
    PhoPz_[nPho_] =  photon.pz();
    PhoPt_[nPho_] =  photon.p4().pt();
    PhoEta_[nPho_] =  photon.p4().eta();
    PhoPhi_[nPho_] =  photon.p4().phi();


    PhoR9_     [nPho_]    =  photon.r9();

/*
    isEBGap[nPho_]    =  photon.isEBGap()? 1:0;
    isEEGap[nPho_]    =  photon.isEEGap()? 1:0;
    isEBEEGap[nPho_]  =  photon.isEBEEGap()? 1:0;
    isTransGap[nPho_] =  (fabs(photon.eta()) > ecalBarrelMaxEta_ && fabs(photon.eta()) < ecalEndcapMinEta_) ? 1:0;
    isEB[nPho_]       =  photon.isEB()? 1:0;
    isEE[nPho_]       =  photon.isEE()? 1:0;
*/

//cout<<"Where_2"<<endl;


// Super-cluster parameters

//    rawEnergy      [nPho_]   =  photon.superCluster()->rawEnergy();
 //   preshowerEnergy[nPho_]   =  photon.superCluster()->preshowerEnergy();
   // numOfPreshClusters[nPho_]=  getNumOfPreshClusters(&photon, e);
    //clustersSize   [nPho_]   =  photon.superCluster()->clustersSize();
    PhoPhiWidth_[nPho_]   =  photon.superCluster()->phiWidth();
    PhoEtaWidth_[nPho_]   =  photon.superCluster()->etaWidth();
    PhoScEta_[nPho_]   =  photon.superCluster()->eta();
    PhoScPhi_[nPho_]   =  photon.superCluster()->phi();
    //scSize         [nPho_]   =  photon.superCluster()->size();

    //ES Ratio
    //PhoESRatio        [nPho_]   =  getESRatio(&photon, e, iSetup);


// Cluster shape variables

    const reco::CaloClusterPtr  seed = photon.superCluster()->seed();

    DetId id     = lazyTool.getMaximum(*seed).first; 
    
    float time  = -999., outOfTimeChi2 = -999., chi2 = -999.;
    int   flags =-1, severity = -1; 
    float seedAppEt = -999;
    float E2E9  = -999.;
    const EcalRecHitCollection & rechits = ( photon.isEB() ? *EBReducedRecHits : *EEReducedRecHits); 
    EcalRecHitCollection::const_iterator it = rechits.find( id );
    if( it != rechits.end() ) { 
	    time = it->time(); 
	    outOfTimeChi2 = it->outOfTimeChi2();
	    chi2 = it->chi2();
	    flags = it->recoFlag();
	    severity = EcalSeverityLevelAlgo::severityLevel( id, rechits, *chStatus );
	    seedAppEt = (id.subdetId() == EcalBarrel)?
	      it->energy()/ cosh( EBDetId::approxEta( id ) ):0;
	    E2E9 = EcalSeverityLevelAlgo::E2overE9( id, rechits, 5.0, 0.0);
    }
/*
    float tlef = -999., tright=-999., ttop=-999., tbottom=-999.;
    std::vector<DetId> left   = lazyTool.matrixDetId(id,-1,-1, 0, 0);
    std::vector<DetId> right  = lazyTool.matrixDetId(id, 1, 1, 0, 0);
    std::vector<DetId> top    = lazyTool.matrixDetId(id, 0, 0, 1, 1);
    std::vector<DetId> bottom = lazyTool.matrixDetId(id, 0, 0,-1,-1);
     cout<<"Where_3"<<endl;

    float *times[4] = {&tlef,&tright,&ttop,&tbottom};
    std::vector<DetId> ids[4]  = {left,right,top,bottom};
    int nt = sizeof(times)/sizeof(float);
    for(int ii=0; ii<nt; ++ii) {
	    if( ids[ii].empty() ) { continue; }
	    it = rechits.find( ids[ii][0] );
	    if( it != rechits.end() ) { *(times[ii]) = it->time(); }
    }
 cout<<"Where_4"<<endl;
  */ 
    PhoSeedTime_[nPho_]  = time;
 
//   seedOutOfTimeChi2     [nPho_]  = outOfTimeChi2;
 //   seedChi2              [nPho_]  = chi2;
 //   seedRecoFlag          [nPho_]  = flags;
 //   seedSeverity          [nPho_]  = severity;
    PhoE2overe9_[nPho_]  = E2E9;
//    seedEt                [nPho_]  = seedAppEt;
    
 //   tLef          [nPho_] = tlef; 
 //   tRight        [nPho_] = tright; 
 //   tTop          [nPho_] = ttop; 
 //   tBottom       [nPho_] = tbottom; 

 //   eMax         [nPho_] =  lazyTool.eMax(*seed);
  //  e2nd         [nPho_] =  lazyTool.e2nd(*seed);
  //  e2x2         [nPho_] =  lazyTool.e2x2(*seed);
  //  e3x2         [nPho_] =  lazyTool.e3x2(*seed);
  //  e3x3         [nPho_] =  lazyTool.e3x3(*seed);
  //  e4x4         [nPho_] =  lazyTool.e4x4(*seed);
  //  e5x5         [nPho_] =  lazyTool.e5x5(*seed);
  //  e2overe8     [nPho_] =  ( lazyTool.e3x3(*seed)-lazyTool.eMax(*seed) ==0 )? 0: lazyTool.e2nd(*seed)/( lazyTool.e3x3(*seed)-lazyTool.eMax(*seed) );

  //  e2x5Right    [nPho_] =  lazyTool.e2x5Right(*seed);
  //  e2x5Left     [nPho_] =  lazyTool.e2x5Left(*seed);
  //  e2x5Top      [nPho_] =  lazyTool.e2x5Top(*seed);
  //  e2x5Bottom   [nPho_] =  lazyTool.e2x5Bottom(*seed);
  //  eRight       [nPho_] =  lazyTool.eRight(*seed);
  //  eLeft        [nPho_] =  lazyTool.eLeft(*seed);
  //  eTop         [nPho_] =  lazyTool.eTop(*seed);
  //  eBottom      [nPho_] =  lazyTool.eBottom(*seed);

  //  vector<float> vCov;
  //  vCov = lazyTool.covariances(*seed);
  //  covPhiPhi    [nPho_] = vCov[2];
  //  covEtaPhi    [nPho_] = vCov[1];
  //  covEtaEta    [nPho_] = vCov[0];

  //  vector<float> viCov;
  //  viCov = lazyTool.localCovariances(*seed);
//cout<<"Here8 ?"<<endl;

    // Photon shower shape parameters 
  //  maxEnergyXtal[nPho_] =  photon.maxEnergyXtal();
 //   sigmaEtaEta  [nPho_] =  photon.sigmaEtaEta();
   PhoSigmaIetaIeta_[nPho_] =  photon.sigmaIetaIeta();
//    sigmaIphiIphi[nPho_] =  sqrt(viCov[2]);
//    sigmaIetaIphi[nPho_] =  sqrt(viCov[1]);
//    r1x5         [nPho_] =  photon.r1x5();
//    r2x5         [nPho_] =  photon.r2x5();
//    e1x5         [nPho_] =  photon.e1x5();
//    e2x5         [nPho_] =  photon.e2x5();



// AOD isolation and identification

    PhohadronicOverEm_[nPho_]   =  photon.hadronicOverEm();
//    hadronicDepth1OverEm[nPho_]   =  photon.hadronicDepth1OverEm();
 //   hadronicDepth2OverEm[nPho_]   =  photon.hadronicDepth2OverEm();
 //   trackIso            [nPho_]   =  photon.trackIso();
 //   caloIso             [nPho_]   =  photon.caloIso();
 //   ecalIso             [nPho_]   =  photon.ecalIso();
 //   hcalIso             [nPho_]   =  photon.hcalIso();

  
    // number of pixel hits                                                                                                                                                             
//    nPixelHits         (nPho_) = nPixel;
    
        
// Delta R= 0.4
    
    PhoecalRecHitSumEtConeDR04_[nPho_]   =  photon.ecalRecHitSumEtConeDR04();
    PhohcalTowerSumEtConeDR04_[nPho_]   =  photon.hcalTowerSumEtConeDR04();
//    hcalDepth1TowerSumEtConeDR04[nPho_]   =  photon.hcalDepth1TowerSumEtConeDR04();
 //   hcalDepth2TowerSumEtConeDR04[nPho_]   =  photon.hcalDepth2TowerSumEtConeDR04();
 //   trkSumPtSolidConeDR04       [nPho_]   =  photon.trkSumPtSolidConeDR04();
    PhotrkSumPtHollowConeDR04_[nPho_]   =  photon.trkSumPtHollowConeDR04();
//    nTrkSolidConeDR04           [nPho_]   =  photon.nTrkSolidConeDR04();
//    nTrkHollowConeDR04          [nPho_]   =  photon.nTrkHollowConeDR04();

   
// Conversion

    //hasConversionTracks[nPho_]   =  photon.hasConversionTracks();
    PhohasPixelSeed_[nPho_]   =  photon.hasPixelSeed();
        

      if (isMC_) {

      edm::Handle<reco::GenParticleCollection> genParticles;

      //  get generated particles and store generator ntuple 
      try { e.getByLabel( genParticlesProducerC,      genParticles );} catch (...) {;}

      float delta(0.15);
      

      const reco::Candidate *cndMc(0);
      reco::GenParticleCollection::const_iterator matchedPart;
      Float_t currentMaxPt(-1);
      
      for (reco::GenParticleCollection::const_iterator it_gen = 
	     genParticles->begin(); it_gen!= genParticles->end(); it_gen++){   

	const reco::Candidate &p = (*it_gen);    
	if (p.status() != 1 || (p.pdgId()) != pdgId_ ) continue;      
	if(ROOT::Math::VectorUtil::DeltaR(p.p4(),phoItr->p4())<delta && p.pt() > currentMaxPt ) {
	  if( p.numberOfMothers() > 0 ) {
	    PhogenMomId_[nPho_] = p.mother()->pdgId();
	    PhogenNSiblings_[nPho_] = p.mother()->numberOfDaughters();
	    if( p.mother()->numberOfMothers() > 0 ) { 
	      PhogenGrandMomId_[nPho_] = p.mother()->mother()->pdgId();
	    }
	  }
	  PhoisGenMatched_[nPho_] = kTRUE; cndMc = &p;
	  currentMaxPt = p.pt();
	  matchedPart  = it_gen;
	}
      }      
	
    	
      // if no matching photon was found try with other particles
      if( ! PhoisGenMatched_[nPho_] ) {

	currentMaxPt = -1;
	for (reco::GenParticleCollection::const_iterator it_gen = 
	       genParticles->begin(); it_gen!= genParticles->end(); it_gen++){
	  const reco::Candidate &p = (*it_gen);    

	  if (p.status() != 1 || find(otherPdgIds_.begin(),otherPdgIds_.end(),fabs(p.pdgId())) == otherPdgIds_.end() ) continue;      	
	  if(ROOT::Math::VectorUtil::DeltaR(p.p4(),phoItr->p4())<delta && p.pt() > currentMaxPt ) {
	    PhogenMomId_[nPho_] = p.pdgId();	
	    if( p.numberOfMothers() > 0 ) {
	      PhogenGrandMomId_[nPho_] = p.mother()->pdgId();
	      PhogenNSiblings_[nPho_]  = p.mother()->numberOfDaughters();
	    }
	    cndMc = &p; // do not set the isGenMatched in this case
	    currentMaxPt = p.pt();
	    matchedPart  = it_gen;
	  }	
	
	} // end of loop over gen particles
      } // if not matched to gen photon

      if(cndMc) {
//	PhogenMatchedP4_ [nPho_]  = TLorentzVector(cndMc->px(),cndMc->py(),cndMc->pz(),cndMc->energy());
	PhogenMatchedPt_[nPho_]  = cndMc->pt();
	PhogenMatchedEta_[nPho_]  = cndMc->eta();
	PhogenMatchedPhi_[nPho_]  = cndMc->phi();
	
//	genIsoDR03_[nPho_]   = getGenCalIso(genParticles,matchedPart,0.3,true,true);
//	genCalIsoDR03_[nPho_]= getGenCalIso(genParticles,matchedPart,0.3,muInGenCaloIso_);
//	genTrkIsoDR03[nPho_]= getGenTrkIso(genParticles,matchedPart,0.3);
	PhogenIsoDR04_[nPho_]   = getGenCalIso(genParticles,matchedPart,0.4,true,true);
	PhogenCalIsoDR04_[nPho_]= getGenCalIso(genParticles,matchedPart,0.4,muInGenCaloIso_);
	PhogenTrkIsoDR04_[nPho_]= getGenTrkIso(genParticles,matchedPart,0.4);
      }

    } // if it's a MC

    //cout<<"Where_3"<<endl;


    
    nPho_++;
   //cout<<"nPho_"<<endl; 
   if (nPho_>kMaxPhotons-1) break;
  }


bool hasPhoReturn = false ;

if(nPho_!=0) hasPhoReturn = true ;
/*
 if(hasPhoReturn)
 {
 cout<<"Here14 ? photon : "<< nPho_ <<" true " << endl;
 }
 else
 {
 cout<<"Here14 ? photon : "<< nPho_ <<" False " << endl;
 }
*/ 
  return (nPho_);
  
}



Float_t PhotonJetAna::getGenCalIso(edm::Handle<reco::GenParticleCollection> handle,
					   reco::GenParticleCollection::const_iterator thisPho,const Float_t dRMax,
					   bool removeMu, bool removeNu)
{
  const Float_t etMin = 0.0;
  Float_t genCalIsoSum = 0.0;
  if(!isMC_)return genCalIsoSum;
  if(!handle.isValid())return genCalIsoSum;

  for (reco::GenParticleCollection::const_iterator it_gen = 
	 handle->begin(); it_gen!=handle->end(); it_gen++){

    if(it_gen == thisPho)continue;      // can't be the original photon
    if(it_gen->status()!=1)continue;    // need to be a stable particle
    if (thisPho->collisionId() != it_gen->collisionId())  // has to come from the same collision
       continue; 
   
    Int_t pdgCode = abs(it_gen->pdgId());
    /// if(pdgCode>11 && pdgCode < 20)continue;     // we should not count neutrinos, muons
    if( removeMu && pdgCode == 13 ) continue;
    if( removeNu && ( pdgCode == 12 || pdgCode == 14 || pdgCode == 16 ) ) continue;

    Float_t et = it_gen->et();
    if(et < etMin) continue; // pass a minimum et threshold, default 0

    Float_t dR = reco::deltaR(thisPho->momentum(), 
			      it_gen->momentum());
    if(dR > dRMax) continue; // within deltaR cone
    genCalIsoSum += et;
    
  }// end of loop over gen particles

  return genCalIsoSum;
}


//=============================================================================
// default cut value of ptMin is 0.0

Float_t PhotonJetAna::getGenTrkIso(edm::Handle<reco::GenParticleCollection> handle,
					   reco::GenParticleCollection::const_iterator thisPho,                                            const Float_t dRMax)
{
  const Float_t ptMin = 0.0;
  Float_t genTrkIsoSum = 0.0;
  if(!isMC_)return genTrkIsoSum;
  if(!handle.isValid())return genTrkIsoSum;

  for (reco::GenParticleCollection::const_iterator it_gen = 
	 handle->begin(); it_gen!=handle->end(); it_gen++){

    if(it_gen == thisPho)continue;      // can't be the original photon
    if(it_gen->status()!=1)continue;    // need to be a stable particle
    if (thisPho->collisionId() != it_gen->collisionId())  // has to come from the same collision
       continue; 
   
   if(it_gen->charge()==0)continue;    // we should not count neutral particles
   
    Float_t pt = it_gen->pt();
    if(pt < ptMin) continue; // pass a minimum pt threshold, default 0

    Float_t dR = reco::deltaR(thisPho->momentum(), 
			      it_gen->momentum());
    if(dR > dRMax) continue; // within deltaR cone
    genTrkIsoSum += pt;
    
  }// end of loop over gen particles

  return genTrkIsoSum;
}









// ------------ method called once each job just before starting event loop  ------------
void 
PhotonJetAna::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonJetAna::endJob() 
{
}
void PhotonJetAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

//define this as a plug-innalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
 }

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonJetAna);
