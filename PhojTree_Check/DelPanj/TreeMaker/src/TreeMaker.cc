// -*- C++ -*-
//
// Package:    TreeMaker
// Class:      TreeMaker
// Original Author:  Anil Singh, Ashok Kumar
//                   Panjab University, Delhi University
//         Created:  Tue Jul  6 21:04:59 CEST 2010
//         Modifid for photon jet  by Yun-Ju Lu 2011 Oct 5        


// system include files
#include <memory>
#include <string>

#include "DelPanj/TreeMaker/interface/TreeMaker.h"

TreeMaker::TreeMaker(const edm::ParameterSet& iConfig)

{
   
   NprocessTree=0;
   fillEventInfo_=0;
   fillGenInfo_=0;
   fillTrigInfo_=0;
   fillMuonInfo_=0;
   fillElecInfo_=0;
   fillJetInfo_=0;
   fillMetInfo_=0;
   fillPhotInfo_=0; 

   fillEventInfo_ = iConfig.getParameter<bool>("fillEventInfo_");
   fillGenInfo_   = iConfig.getParameter<bool>("fillGenInfo_");
   fillMuonInfo_  = iConfig.getParameter<bool>("fillMuonInfo_");
   fillElecInfo_  = iConfig.getParameter<bool>("fillElecInfo_");
   fillJetInfo_   = iConfig.getParameter<bool>("fillJetInfo_");
   fillMetInfo_   = iConfig.getParameter<bool>("fillMetInfo_");
   fillTrigInfo_  = iConfig.getParameter<bool>("fillTrigInfo_");
   fillPhotInfo_  = iConfig.getParameter<bool>("fillPhotInfo_");


   outFileName_= iConfig.getParameter<std::string>("outFileName");
   file = new TFile(outFileName_.c_str(),"recreate");
;
   tree_ = new TTree("tree","tree");
   if( fillEventInfo_) eventInfo_=new eventInfo("eventInfo",tree_,iConfig); 
   if( fillGenInfo_)   genInfoTree_ = new genInfoTree("gen",tree_,iConfig);
   if( fillMuonInfo_)  patMuTree_= new patMuonTree("mu",tree_,iConfig);
   if( fillElecInfo_)  patElecTree_ = new patElecTree("ele",tree_,iConfig);
   if( fillMetInfo_)   patMetTree_= new patMetTree("met",tree_,iConfig);
   if( fillJetInfo_)   jetTree_=new jetTree("patJetPfAk05",tree_,iConfig);
   if( fillTrigInfo_)  patHltTree_ = new patHltTree("Hlt",tree_,iConfig); 
   if( fillPhotInfo_)  photonTree_ = new photonTree("Photon", tree_, iConfig); 
}


TreeMaker::~TreeMaker()
{
  //delete tree_;
}

void
TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  NprocessTree++;
  using namespace edm;
  //check is Data or MC   
 
  if( fillEventInfo_) eventInfo_   ->Fill(iEvent);
  if( fillGenInfo_)   genInfoTree_ ->Fill(iEvent);
  if( fillElecInfo_)  patElecTree_ ->Fill(iEvent);
  if( fillMuonInfo_)  patMuTree_   ->Fill(iEvent);
  if( fillJetInfo_)   jetTree_     ->Fill(iEvent,iSetup);
  if( fillMetInfo_)   patMetTree_  ->Fill(iEvent);
  if( fillTrigInfo_)  patHltTree_  ->Fill(iEvent,iSetup,hltConfig_);
  if( fillPhotInfo_)  photonTree_  ->Fill(iEvent,iSetup);
  tree_->Fill();
}

void
TreeMaker::beginJob(){
}

void TreeMaker::beginRun(const edm::Run & iRun, const edm::EventSetup & iSetup) {
   bool changed(true);
   if (hltConfig_.init(iRun, iSetup, "HLT", changed)) {
     // if init returns TRUE, initialisation has succeeded!
     if (changed) {
       // The HLT config has actually changed wrt the previous Run, hence rebook your
       // histograms or do anything else dependent on the revised HLT config
     }
   } else {
     // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
     // with the file and/or code and needs to be investigated!
     std::cout << " HLT config extraction failure with process name " << "HLT" << std::endl;
     //LogError("MyAnalyzer") << " HLT config extraction failure with process name " << processName_;
     // In this case, all access methods will return empty values!
   }
 }


void
TreeMaker::endJob() {
  std::cout<<"TreeMaker processes "<<NprocessTree <<std::endl;
  file->cd();
  file->Write();
  file->Close();
}

DEFINE_FWK_MODULE(TreeMaker);
