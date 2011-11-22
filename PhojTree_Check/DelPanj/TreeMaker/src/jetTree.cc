#include "DelPanj/TreeMaker/interface/jetTree.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "TMath.h"
typedef math::XYZTLorentzVector LorentzVector;

jetTree::jetTree(std::string desc, TTree* tree, const edm::ParameterSet& iConfig):
baseTree(desc, tree),
parSet_(iConfig.getParameter<edm::ParameterSet>(desc.c_str()))
{
  JetLabel_ = parSet_.getParameter<edm::InputTag>("Jets");
  SetBranches();
}


jetTree::~jetTree(){
  delete tree_;
}


void
jetTree::Fill(const edm::Event& iEvent){
  Clear();
  
  edm::Handle<std::vector<pat::Jet> > JetHandle;
  if(not iEvent.getByLabel(JetLabel_,JetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
         <<JetLabel_<<std::endl; exit(0);}

  const std::vector<pat::Jet>* jets = JetHandle.product();
  std::vector<pat::Jet>::const_iterator jet =jets->begin();

  for(;jet!=jets->end();jet++){
    //add only save jet with pt >18 eta>3.2 
    if(jet->pt()<18) continue ;
    if(jet->eta()>3.2) continue;

    //Stuff common for all jets.
    JetPt_.push_back(jet->pt());
    JetEta_.push_back(jet->eta());
    JetPhi_.push_back(jet->phi());
    JetM_.push_back(jet->mass());
    //JetRapidity_.push_back(jet->Rap());
    JetPx_.push_back(jet->px());
    JetPy_.push_back(jet->py());
    JetPz_.push_back(jet->pz());
    JetEn_.push_back(jet->energy());
 
    double pt = -99;
    double px = -99;
    double py = -99;
    double pz = -99;
    double en = -99;
    //Check if uncorrected jet information
    //is available.
    /*
     The PAT allows us to know the names for the
     JEC correction that have been applied. The
     code here is only for the debugging and 
     studying the JEC. And Not included in general
     runs of our code.
    */
    std::string jecset = jet->currentJECSet();
    //std::cout<<"JEC Set is: "<<jecset<<std::endl;			 
    bool isAvailable = jet->jecSetAvailable("Uncorrected");
     if(isAvailable){
      LorentzVector uncorrP4 = jet->correctedP4("Uncorrected");
      pt = uncorrP4.pt();
      px = uncorrP4.px();
      py = uncorrP4.py();
      pz = uncorrP4.pz();
      en = uncorrP4.e();
    }
    JetUnCorrPt_.push_back(pt);
    JetUnCorrPx_.push_back(px);
    JetUnCorrPy_.push_back(py);
    JetUnCorrPz_.push_back(pz);
    JetUnCorrEnt_.push_back(en);

    if(jet->isPFJet()){
      
      JetPhotEn_.push_back(jet->photonEnergy());
      JetElecEn_.push_back(jet->electronEnergy());
      JetMuonEn_.push_back(jet->muonEnergy());
      JetHfHadEn_.push_back(jet->HFHadronEnergy());
      JetHfEmEn_.push_back(jet->HFEMEnergy());
      JetCharHadE_.push_back(jet->chargedHadronEnergy ());
      JetNeutHadE_.push_back(jet->neutralHadronEnergy ());   
      JetCharEmE_.push_back(jet->chargedEmEnergy ());
      JetCharMuE_.push_back(jet->chargedMuEnergy ());
      JetNeutEmE_.push_back(jet->neutralEmEnergy ());
      
      JetMuonMulti_.push_back(jet->muonMultiplicity());
      JetNeutMulti_.push_back(jet->neutralMultiplicity ());   
      JetCharMulti_.push_back(jet->chargedMultiplicity());
      JetCharHadMulti_.push_back(jet->chargedHadronMultiplicity());
      JetNeutHadMulti_.push_back(jet->neutralHadronMultiplicity());
      JetPhotMulti_.push_back(jet->photonMultiplicity());
      JetElecMulti_.push_back(jet->electronMultiplicity());
      JetHfHadMulti_.push_back(jet->HFHadronMultiplicity());
      JetHfEmMulti_.push_back(jet->HFEMMultiplicity());
     
      JetPhotEnFr_.push_back(jet->photonEnergyFraction());
      JetMuonEnFr_.push_back(jet->muonEnergyFraction());
      JetHfHadEnFr_.push_back(jet->HFHadronEnergyFraction());
      JetHfEmEnFr_.push_back(jet->HFEMEnergyFraction());
      JetNeutEmEFr_.push_back(jet->neutralEmEnergyFraction ());
      JetCharHadEFr_.push_back(jet->chargedHadronEnergyFraction ());
      JetNeutHadEFr_.push_back(jet->neutralHadronEnergyFraction ());
      JetCharEmEFr_.push_back(jet->chargedEmEnergyFraction ());
      JetCharMuEFr_.push_back(jet->chargedMuEnergyFraction ());
      JetNConstituents_.push_back(jet->getPFConstituents().size()); 
   }
    
    if(jet->isCaloJet()){
      JetN60_.push_back(jet->n60());
      JetN90_.push_back(jet->n90());
      JetEmFr_.push_back(jet->emEnergyFraction());
      JetHadFr_.push_back(jet->energyFractionHadronic());
      JetEmEbEn_.push_back(jet->emEnergyInEB());
      JetEmEeEn_.push_back(jet->emEnergyInEE());
      JetEmHfEn_.push_back(jet->emEnergyInHF());
      JetHadHbEn_.push_back(jet->hadEnergyInHB());
      JetHadHeEn_.push_back(jet->hadEnergyInHE());
      JetHadHfEn_.push_back(jet->hadEnergyInHF());
      JetHadHoEn_.push_back(jet->hadEnergyInHO());    
    }

  }
}

void
jetTree::SetBranches(){
  
  AddBranch(&JetPt_, "Pt_");
  AddBranch(&JetEta_, "Eta_");
  AddBranch(&JetPhi_, "Phi_");
  AddBranch(&JetM_, "M_");
  AddBranch(&JetRapidity_, "Rapidity_");
  AddBranch(&JetPx_, "Px_");
  AddBranch(&JetPy_, "Py_");
  AddBranch(&JetPz_, "Pz_");
  AddBranch(&JetEn_, "En_");
  AddBranch(&JetUnCorrPt_, "UnCorrPt_");
  AddBranch(&JetUnCorrPx_, "UnCorrPx_");
  AddBranch(&JetUnCorrPy_, "UnCorrPy_");
  AddBranch(&JetUnCorrPz_, "UnCorrPz_");
  AddBranch(&JetUnCorrEnt_, "UnCorrEnt_");
  
  //Add PFlow Specific Branches.
  //if(storePFlow_){
    AddBranch(&JetPhotEn_, "PhotEn_");
    AddBranch(&JetElecEn_, "ElecEn_");
    AddBranch(&JetMuonEn_, "MuonEn_");
    AddBranch(&JetHfHadEn_, "HfHadEn_");
    AddBranch(&JetHfEmEn_, "HfEmEn_");
    AddBranch(&JetCharHadE_, "CharHadE_");
    AddBranch(&JetNeutHadE_, "NeutHadE_");
    AddBranch(&JetCharEmE_, "CharEmE_");
    AddBranch(&JetCharMuE_, "CharMuE_");
    AddBranch(&JetNeutEmE_, "NeutEmE_");
    AddBranch(&JetMuonMulti_, "MuonMulti_");
    AddBranch(&JetNeutMulti_, "NeutMulti_");
    AddBranch(&JetCharMulti_, "CharMulti_");
    AddBranch(&JetCharHadMulti_, "CharHadMulti_");
    AddBranch(&JetNeutHadMulti_, "NeutHadMulti_");
    AddBranch(&JetPhotMulti_, "PhotMulti_");
    AddBranch(&JetElecMulti_, "ElecMulti_");
    AddBranch(&JetHfHadMulti_, "HfHadMulti_");
    AddBranch(&JetHfEmMulti_, "HfEmMulti_");
    AddBranch(&JetPhotEnFr_, "PhotEnFr_");
    AddBranch(&JetMuonEnFr_, "MuonEnFr_");
    AddBranch(&JetHfHadEnFr_, "HfHadEnFr_");
    AddBranch(&JetHfEmEnFr_, "HfEmEnFr_");
    AddBranch(&JetNeutEmEFr_, "NeutEmEFr_");
    AddBranch(&JetCharHadEFr_, "CharHadEFr_");
    AddBranch(&JetNeutHadEFr_, "NeutHadEFr_");
    AddBranch(&JetCharEmEFr_, "CharEmEFr_");
    AddBranch(&JetCharMuEFr_, "CharMuEFr_");
    AddBranch(&JetNConstituents_, "JetNConstituents_");
  //}
  //Add CaloSpecific Branches
  //if(storeCalo_){
/*
    AddBranch(&JetN60_,      "N60_");
    AddBranch(&JetN90_,      "N90_");
    AddBranch(&JetEmFr_,     "EmFr_");
    AddBranch(&JetHadFr_,    "HadFr_");
    AddBranch(&JetEmEbEn_,   "EmEbEn_");
    AddBranch(&JetEmEeEn_,   "EmEeEn_");
    AddBranch(&JetEmHfEn_,   "EmHfEn_");
    AddBranch(&JetHadHbEn_,  "HadHbEn_");
    AddBranch(&JetHadHeEn_,  "HadHeEn_");
    AddBranch(&JetHadHfEn_,  "HadHfEn_");
    AddBranch(&JetHadHoEn_,  "HadHoEn_");
    AddBranch(&JetTotalEm_,  "TotalEm_");
    AddBranch(&JetTotalHad_, "TotalHad_");
    AddBranch(&JetEmEbFr_,   "EmEbFr_");
    AddBranch(&JetEmEeFr_,   "EmEeFr_");
    AddBranch(&JetEmHfFr_,   "EmHfFr_");
    AddBranch(&JetHadHbFr_,  "HadHbFr_");
    AddBranch(&JetHadHeFr_,  "HadHeFr_");
    AddBranch(&JetHadHfFr_,  "HadHfFr_");
    AddBranch(&JetHadHoFr_,  "HadHoFr_");
    //AddBranch(&JetHadHoFr_,"HadHoFr_");
 */ 
 //}
}


void
jetTree::Clear(){
  JetPt_.clear();
  JetEta_.clear();
  JetPhi_.clear();
  JetM_.clear();
  JetRapidity_.clear();
  JetPx_.clear();
  JetPy_.clear();
  JetPz_.clear();
  JetEn_.clear();
  JetUnCorrPt_.clear();
  JetUnCorrPx_.clear();
  JetUnCorrPy_.clear();
  JetUnCorrPz_.clear();
  JetUnCorrEnt_.clear();
  JetPhotEn_.clear();
  JetElecEn_.clear();
  JetMuonEn_.clear();
  JetHfHadEn_.clear();
  JetHfEmEn_.clear();
  JetCharHadE_.clear();
  JetNeutHadE_.clear();
  JetCharEmE_.clear();
  JetCharMuE_.clear();
  JetNeutEmE_.clear();
  JetMuonMulti_.clear();
  JetNeutMulti_.clear();
  JetCharMulti_.clear();
  JetCharHadMulti_.clear();
  JetNeutHadMulti_.clear();
  JetPhotMulti_.clear();
  JetElecMulti_.clear();
  JetHfHadMulti_.clear();
  JetHfEmMulti_.clear();
  JetPhotEnFr_.clear();
  JetMuonEnFr_.clear();
  JetHfHadEnFr_.clear();
  JetHfEmEnFr_.clear();
  JetNeutEmEFr_.clear();
  JetCharHadEFr_.clear();
  JetNeutHadEFr_.clear();
  JetCharEmEFr_.clear();
  JetCharMuEFr_.clear();
  JetN60_.clear();
  JetN90_.clear();
  JetEmFr_.clear();
  JetHadFr_.clear();
  JetEmEbEn_.clear();
  JetEmEeEn_.clear();
  JetEmHfEn_.clear();
  JetHadHbEn_.clear();
  JetHadHeEn_.clear();
  JetHadHfEn_.clear();
  JetHadHoEn_.clear();
  JetTotalEm_.clear(); 
  JetTotalHad_.clear();
  JetEmEbFr_.clear();
  JetEmEeFr_.clear();
  JetEmHfFr_.clear();
  JetHadHbFr_.clear();
  JetHadHeFr_.clear();
  JetHadHfFr_.clear();
  JetHadHoFr_.clear();
  JetNConstituents_.clear();
}
