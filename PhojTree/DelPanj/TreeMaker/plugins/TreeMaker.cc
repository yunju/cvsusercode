// -*- C++ -*-
//
// Package:    TreeMaker
// Class:      TreeMaker
// Original Author:  Anil Singh, Ashok Kumar
//                   Panjab University, Delhi University
//         Created:  Tue Jul  6 21:04:59 CEST 2010
//


// system include files
#include <memory>
#include <string>

#include "DelPanj/TreeMaker/interface/TreeMaker.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

TreeMaker::TreeMaker(const edm::ParameterSet& iConfig)

{

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


   //outFileName_= iConfig.getParameter<std::string>("outFileName");
   edm::Service<TFileService> fs;

   //file = new TFile(outFileName_.c_str(),"recreate");
;
   //tree_ = new TTree("tree","tree");

   tree_ = fs->make<TTree>("tree","tree");
   if( fillEventInfo_) eventInfo_=new eventInfo("eventInfo",tree_); 
   if( fillGenInfo_)   genInfoTree_ = new genInfoTree("gen",tree_,iConfig);
   if( fillMuonInfo_)  patMuTree_= new patMuonTree("mu",tree_,iConfig);
   if( fillElecInfo_)  patElecTree_ = new patElecTree("ele",tree_,iConfig);
   if( fillMetInfo_)   patMetTree_= new patMetTree("met",tree_,iConfig);
   if( fillJetInfo_)   jetTree_=new jetTree("patJetPfAk05",tree_,iConfig);
   if( fillTrigInfo_)  patHltTree_ = new patHltTree("Hlt",tree_); 
   if( fillPhotInfo_)  photonTree_ = new photonTree("patPhoton", tree_, iConfig); 
}


TreeMaker::~TreeMaker()
{
  //delete tree_;
}

void
TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  
  if( fillEventInfo_) eventInfo_   ->Fill(iEvent);
  if( fillGenInfo_)   genInfoTree_ ->Fill(iEvent);
  if( fillElecInfo_)  patElecTree_ ->Fill(iEvent);
  if( fillMuonInfo_)  patMuTree_   ->Fill(iEvent);
  if( fillJetInfo_)   jetTree_     ->Fill(iEvent);
  if( fillMetInfo_)   patMetTree_  ->Fill(iEvent);
  if( fillTrigInfo_)  patHltTree_  ->Fill(iEvent);
  if( fillPhotInfo_)  photonTree_  ->Fill(iEvent);
  tree_->Fill();
}

void
TreeMaker::beginJob(){
}

void
TreeMaker::endJob() {
//  file->cd();
//  file->Write();
//  file->Close();
}

DEFINE_FWK_MODULE(TreeMaker);
