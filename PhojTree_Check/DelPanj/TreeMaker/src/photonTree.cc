//#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
//#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DelPanj/TreeMaker/interface/photonTree.h"

photonTree::photonTree(std::string name, TTree* tree, const pset& iConfig):baseTree(name,tree){
  photonLabel_  = iConfig.getParameter<edm::InputTag> ("photonLabel");
  reducedEcalRecHitsEBLabel_  = iConfig.getParameter<edm::InputTag> ("reducedEcalRecHitsEBLabel");
  reducedEcalRecHitsEELabel_  = iConfig.getParameter<edm::InputTag> ("reducedEcalRecHitsEELabel");
  genParticlesProducerC       = iConfig.getParameter<edm::InputTag> ("genParticlesProducerPY");
  usePFObjects_ = iConfig.getParameter<bool> ("usePFlow");
  pdgId_        = iConfig.getUntrackedParameter<int>("pdgId", 22);   
  otherPdgIds_                     = iConfig.getUntrackedParameter<vector<int> >("OtherPdgIds", vector<int>(1,11) );
  muInGenCaloIso_                  = iConfig.getUntrackedParameter<int>("muInGenCaloIso", 10); 
  isMC_   = iConfig.getParameter<bool>("isMCPY");
  rho25C           = iConfig.getParameter<InputTag>("rho25PY");  
  rho44C           = iConfig.getParameter<InputTag>("rho44PY");  
 


  SetBranches();
}

photonTree::~photonTree(){
  delete tree_;
}

void photonTree::Fill(const edm::Event& iEvent,const edm::EventSetup& iSetup){
  Clear();
  
  //check if input is data or MC
   
   
 //  std::cout<<"isMC_"<<isMC_<<std::endl;
  
  //saveing fast jet rho correction
    edm::Handle<double> rhoHandle25;
    iEvent.getByLabel(rho25C, rhoHandle25);
 
    rho25_ = *(rhoHandle25.product()); 

    edm::Handle<double> rhoHandle44;
    iEvent.getByLabel(rho44C, rhoHandle44);
   
    rho44_ = *(rhoHandle44.product());    


  //fetch the input collection
  edm::Handle<pat::PhotonCollection> photonHandle;
  if(not iEvent.getByLabel(photonLabel_,photonHandle)){
    std::cout<<"FATAL EXCEPTION(photon tree): "<<"Following Not Found: "<<photonLabel_<<std::endl; 
    exit(0);
  }  
  pat::PhotonCollection phColl(*(photonHandle.product()));

//sort the objects by transverse momentum
 
  GreaterByPt<pat::Photon> pTComparator_;
  std::sort(phColl.begin(), phColl.end(), pTComparator_);

  nPhoton_ = phColl.size();
 //using cluster lazytool
 
  edm::Handle<EcalRecHitCollection> EBReducedRecHits;
  iEvent.getByLabel(reducedEcalRecHitsEBLabel_, EBReducedRecHits);
  edm::Handle<EcalRecHitCollection> EEReducedRecHits;
  iEvent.getByLabel(reducedEcalRecHitsEELabel_, EEReducedRecHits);

 

  EcalClusterLazyTools lazyTool(iEvent, iSetup,reducedEcalRecHitsEBLabel_,reducedEcalRecHitsEELabel_);   
  

  pat::PhotonCollection::const_iterator ph;
  
  for(ph=phColl.begin(); ph!=phColl.end(); ph++){
  pat::Photon photon = pat::Photon(*ph);
    photonPt_.push_back(photon.pt());
    photonEta_.push_back(photon.p4().eta());
    photonPhi_.push_back(photon.p4().phi());
    photonEt_.push_back(photon.et()); 
    photonEnergy_.push_back(photon.energy());
    photonPx_.push_back(photon.px());
    photonPy_.push_back(photon.py());
    photonPz_.push_back(photon.pz());
    photonR9_.push_back(photon.r9());
    photonPhiWidth_.push_back(photon.superCluster()->phiWidth());
    photonEtaWidth_.push_back(photon.superCluster()->etaWidth());

    photonScPhi_.push_back(photon.superCluster()->phi());

    photonScEta_.push_back(photon.superCluster()->eta());
    
    
    photonSigmaIetaIeta_.push_back(photon.sigmaIetaIeta());
    //if(photon.pt()>75&&fabs(photon.p4().eta())<2.5&&photon.hadronicOverEm()<0.05)
   // cout<<"Tree test "<<photon.pt() <<" "<<photon.p4().eta() <<photon.sigmaIetaIeta()<<" "<< photon.ecalRecHitSumEtConeDR04()<<" "<<photon.hcalTowerSumEtConeDR04()<<" "<<photon.trkSumPtHollowConeDR04()<<" "<<photon.hasPixelSeed()<<" "<<photon.hadronicOverEm()<<endl;
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
	   // severity = EcalSeverityLevelAlgo::severityLevel( id, rechits, *chStatus );
	    seedAppEt = (id.subdetId() == EcalBarrel)?
	      it->energy()/ cosh( EBDetId::approxEta( id ) ):0;
	    //E2E9 = EcalSeverityLevelAlgo::E2overE9( id, rechits, 5.0, 0.0);
    }




    photonSigmaIphiIphi_.push_back(-999);
   
    photonSeedTime_.push_back(time);
    photonseedSeverity_.push_back(-999);
    photonE2overe9_.push_back(-999);
   
    photonhadronicOverEm_.push_back(photon.hadronicOverEm());
    photonecalRecHitSumEtConeDR04_.push_back(photon.ecalRecHitSumEtConeDR04());
    photonhcalTowerSumEtConeDR04_.push_back(photon.hcalTowerSumEtConeDR04());
    photontrkSumPtHollowConeDR04_.push_back(photon.trkSumPtHollowConeDR04());
    photonhasPixelSeed_.push_back(photon.hasPixelSeed());
    photonESRatio_.push_back(ph->et());
   
    

    if (isMC_) {

      edm::Handle<reco::GenParticleCollection> genParticles;
      
      if(not iEvent.getByLabel(genParticlesProducerC,genParticles)){
        std::cout<<"FATAL EXCEPTION(photon tree gen): "<<"Following Not Found: "<<genParticlesProducerC<<std::endl; 
        exit(0);
      }  
      edm::Handle<GenEventInfoProduct>    genEventScale;
      if(not iEvent.getByLabel("generator", genEventScale)){
      std::cout<<"FATAL EXCEPTION(photon tree generator): "<<"Following Not Found: "<<"generator"<<std::endl;
      exit(0);
      }
      MCpthat_ = genEventScale->binningValues()[0];
    

      float delta(0.15);
      

      const reco::Candidate *cndMc(0);
      reco::GenParticleCollection::const_iterator matchedPart;
      Float_t currentMaxPt(-1);
     // cout<<PhoisGenMatched_[nPho_]<<endl; 
      bool PhoGenMatched = false ; // intiialize GenMatch
      int  ThePhogenMomId= -9999;
      int  ThePhogenNSiblings =-9999;
      int  ThegenGrandMomId =-9999;  
      //  cout<<PhoisGenMatched_[nPho_]<<endl; 
      for (reco::GenParticleCollection::const_iterator it_gen = 
	     genParticles->begin(); it_gen!= genParticles->end(); it_gen++){   
      

	const reco::Candidate &p = (*it_gen);    
	if (p.status() != 1 || (p.pdgId()) != pdgId_ ) continue;      
	if(ROOT::Math::VectorUtil::DeltaR(p.p4(),ph->p4())<delta && p.pt() > currentMaxPt ) {
	   if( p.numberOfMothers() == 1 ) {
	    ThePhogenMomId=p.mother()->pdgId();
            ThePhogenNSiblings =p.mother()->numberOfDaughters();
	    if( p.mother()->numberOfMothers() ==1 ) { 
              ThegenGrandMomId =p.mother()->mother()->pdgId();
	    }
	  }
	  PhoGenMatched  = true; cndMc = &p;
	  currentMaxPt = p.pt();
	  matchedPart  = it_gen;
	}
      }//gen loop      
	
      // if no matching photon was found try with other particles
      if( ! PhoGenMatched ) {

	currentMaxPt = -1;
	for (reco::GenParticleCollection::const_iterator it_gen = 
	       genParticles->begin(); it_gen!= genParticles->end(); it_gen++){
	  const reco::Candidate &p = (*it_gen);    

	  if (p.status() != 1 || find(otherPdgIds_.begin(),otherPdgIds_.end(),fabs(p.pdgId())) == otherPdgIds_.end() ) continue;      	
	  if(ROOT::Math::VectorUtil::DeltaR(p.p4(),ph->p4())<delta && p.pt() > currentMaxPt ) {
 
	     if( p.numberOfMothers() ==1 ) {

              ThePhogenMomId=p.mother()->pdgId();
              ThePhogenNSiblings =p.mother()->numberOfDaughters();
	    }
	    cndMc = &p; // do not set the isGenMatched in this case
            currentMaxPt = p.pt();
	    matchedPart  = it_gen;
	  }	
	
	} // end of loop over gen particles
      } // if not matched to gen photon
      PhoisGenMatched_.push_back(PhoGenMatched);	
      PhogenMomId_.push_back(ThePhogenMomId); 
      PhogenNSiblings_.push_back(ThePhogenNSiblings) ;
      PhogenGrandMomId_.push_back(ThegenGrandMomId);
     
      if(cndMc) {
//	PhogenMatchedP4_ [nPho_]  = TLorentzVector(cndMc->px(),cndMc->py(),cndMc->pz(),cndMc->energy());
	PhogenMatchedE_ .push_back(cndMc->energy());
        PhogenMatchedPx_.push_back(cndMc->px());
        PhogenMatchedPy_ .push_back(cndMc->py());
        PhogenMatchedPz_.push_back(cndMc->pz());
        PhogenMatchedPt_.push_back(cndMc->pt());
	PhogenMatchedEta_.push_back(cndMc->eta());
	PhogenMatchedPhi_.push_back(cndMc->phi());
	
//	genIsoDR03_[nPho_]   = getGenCalIso(genParticles,matchedPart,0.3,true,true);
//	genCalIsoDR03_[nPho_]= getGenCalIso(genParticles,matchedPart,0.3,muInGenCaloIso_);
//	genTrkIsoDR03[nPho_]= getGenTrkIso(genParticles,matchedPart,0.3);
	PhogenIsoDR04_ .push_back(getGenCalIso(genParticles,matchedPart,0.4,true,true));
	PhogenCalIsoDR04_.push_back( getGenCalIso(genParticles,matchedPart,0.4,muInGenCaloIso_));
	PhogenTrkIsoDR04_.push_back( getGenTrkIso(genParticles,matchedPart,0.4));
      }






    } // if it's a MC
  }  

}

void photonTree::SetBranches(){
  AddBranch(&rho25_  ,"rho25");
  AddBranch(&rho44_ ,"rho44");
  AddBranch(&nPhoton_  ,"NumPh_");
  AddBranch(&photonPt_ ,"Pt");
  AddBranch(&photonEta_,"Eta");
  AddBranch(&photonPhi_,"Phi");
  AddBranch(&photonEt_ , "Et");
  AddBranch(&photonEnergy_,"Energy");
  AddBranch(&photonPx_,"Px");
  AddBranch(&photonPy_,"Py");
  AddBranch(&photonPz_,"Pz");
  AddBranch(&photonR9_,"R9");
  AddBranch(&photonPhiWidth_,"PhiWidth");
  AddBranch(&photonEtaWidth_,"EtaWidth");
  AddBranch(&photonScPhi_,"ScPhi");
  AddBranch(&photonScEta_,"ScEta");
 // AddBranch(&photonESRatio_,"ESRatio");
  AddBranch(&photonSigmaIetaIeta_,"SigmaIetaIeta");
 // AddBranch(&photonSigmaIphiIphi_,"SigmaIphiIphi");
  AddBranch(&photonSeedTime_,"SeedTime");
  AddBranch(&photonseedSeverity_,"seedSeverity");
 // AddBranch(&photonE2overe9_,"E2overe9");
  AddBranch(&photonhadronicOverEm_,"hadronicOverEm");
  AddBranch(&photonecalRecHitSumEtConeDR04_,"ecalRecHitSumEtConeDR04");
  AddBranch(&photonhcalTowerSumEtConeDR04_,"hcalTowerSumEtConeDR04");
  AddBranch(&photontrkSumPtHollowConeDR04_,"trkSumPtHollowConeDR04");
  AddBranch(&photonhasPixelSeed_,"hasPixelSeed");

//  std::cout<<"isMC_Set"<<isMC_<<std::endl;
//cin.get();


  if (isMC_) {

     AddBranch(&MCpthat_,"MCpthat");
     AddBranch(&PhoisGenMatched_,"isGenMatched");
     AddBranch(&PhogenMomId_,"genMomId"); 
     AddBranch(&PhogenGrandMomId_,"genGrandMomId");
     AddBranch(&PhogenNSiblings_,"genNSiblings");
     AddBranch(&PhogenMatchedE_,"genMatchedE");
     AddBranch(&PhogenMatchedPx_,"genMatchedPx");
     AddBranch(&PhogenMatchedPy_,"genMatchedPy");
     AddBranch(&PhogenMatchedPz_,"genMatchedPz");
     AddBranch(&PhogenMatchedPt_,"genMatchedPt");
     AddBranch(&PhogenMatchedEta_,"genMatchedEta");
     AddBranch(&PhogenMatchedPhi_,"genMatchedPhi");
     AddBranch(&PhogenCalIsoDR04_,"genCalIsoDR04");
     AddBranch(&PhogenTrkIsoDR04_,"genTrkIsoDR04");
     AddBranch(&PhogenIsoDR04_,"genIsoDR04");
}

}

void photonTree::Clear(){
  photonPt_.clear();
  photonEta_.clear();
  photonPhi_.clear();
  photonEt_.clear();
  photonEnergy_.clear();
  photonPx_.clear();
  photonPy_.clear();
  photonPz_.clear();
  photonR9_.clear();
  photonPhiWidth_.clear();
  photonEtaWidth_.clear();
  photonScPhi_.clear();
  photonScEta_.clear();
  photonESRatio_.clear();
  photonSigmaIetaIeta_.clear();
  photonSigmaIphiIphi_.clear();
  photonSeedTime_.clear();
  photonseedSeverity_.clear();
  photonE2overe9_.clear();
  photonhadronicOverEm_.clear();
  photonecalRecHitSumEtConeDR04_.clear();
  photonhcalTowerSumEtConeDR04_.clear();
  photontrkSumPtHollowConeDR04_.clear();
  photonhasPixelSeed_.clear();
  PhoisGenMatched_.clear();
  PhogenMomId_.clear(); 
  PhogenGrandMomId_.clear();
  PhogenNSiblings_.clear();
  PhogenMatchedE_.clear();
  PhogenMatchedPx_.clear();
  PhogenMatchedPy_.clear();
  PhogenMatchedPz_.clear();
  PhogenMatchedPt_.clear();
  PhogenMatchedEta_.clear();
  PhogenMatchedPhi_.clear();
  PhogenCalIsoDR04_.clear();
  PhogenTrkIsoDR04_.clear();
  PhogenIsoDR04_.clear();


   







  }
Float_t photonTree::getGenCalIso(edm::Handle<reco::GenParticleCollection> handle,
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

Float_t photonTree::getGenTrkIso(edm::Handle<reco::GenParticleCollection> handle,
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





