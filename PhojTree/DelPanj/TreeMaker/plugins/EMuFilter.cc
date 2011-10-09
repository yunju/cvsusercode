// -*- C++ -*-
//
// Package:    EMuFilter
// Class:      EMuFilter
// 
/**\class EMuFilter EMuFilter.cc DelPanj/EMuFilter/src/EMuFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lovedeep Kaur (Panjab U)
//         Created:  Mon Sep 19 14:33:57 CDT 2011
// $Id: EMuFilter.cc,v 1.1 2011/10/05 22:15:29 anil Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/eSelector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "DelPanj/TreeMaker/interface/muSelector.h"
#include "DelPanj/TreeMaker/interface/muHist.h"
#include "DelPanj/TreeMaker/interface/eHist.h"

//
//
// class declaration
//

class EMuFilter : public edm::EDFilter {
   public:
      explicit EMuFilter(const edm::ParameterSet&);
      ~EMuFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

    //electron tag
     edm::InputTag patElecLabel_;
     edm::InputTag patMuLabel_;
     edm::InputTag beamSpotLabel_;


  //electron selection cuts.
     eSelector ewp_;//for lead ele.
     muSelector mwp_;//for sublead ele
     double zMassLow_;
     double zMassHigh_;


     bool scoreHistos_;	
     eHist* e1Hist_;
     muHist* m1Hist_;
     

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EMuFilter::EMuFilter(const edm::ParameterSet& iConfig):
patMuLabel_(iConfig.getParameter<edm::InputTag>("patMuons")),
beamSpotLabel_(iConfig.getParameter<edm::InputTag>("beamSpot")),
patElecLabel_(iConfig.getParameter<edm::InputTag>("patElectrons")),
ewp_(iConfig.getParameter<edm::ParameterSet>("elPset")),
mwp_(iConfig.getParameter<edm::ParameterSet>("muPset")),
zMassLow_(iConfig.getParameter<double>("zMassLowerLimit")),
zMassHigh_(iConfig.getParameter<double>("zMassUpperLimit")),
scoreHistos_(iConfig.getParameter<double>("scoreHistos"))
{
   //now do what ever initialization is needed
    edm::Service<TFileService> fs;
    if(scoreHistos_){
      e1Hist_ = new eHist("LeadingElectron", fs);
      m1Hist_ = new muHist("leadingMuon", fs);
     }

}


EMuFilter::~EMuFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
EMuFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{ //std::cout<<"test2\n";
 using namespace edm;
   edm::Handle<pat::ElectronCollection> patElecHandle;
     if(not iEvent.getByLabel(patElecLabel_,patElecHandle)){
       std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
                <<patElecLabel_<<std::endl; return false;
     }
     pat::ElectronCollection eleColl(*(patElecHandle.product()));
     if(eleColl.size()<1)return false;


     edm::Handle<pat::MuonCollection> patMuHandle;
     if(not iEvent.getByLabel(patMuLabel_,patMuHandle)){
       std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
        <<patMuLabel_<<std::endl; return false;
     }
     pat::MuonCollection muColl(*(patMuHandle.product()));
     if(muColl.size()<1)return false;
    

     edm::Handle<reco::BeamSpot> beamSpotHandle;
     iEvent.getByLabel(beamSpotLabel_, beamSpotHandle);
     reco::BeamSpot beamSpot = (*beamSpotHandle);
 	

     pat::Electron l1 = eleColl[0];
     pat::Muon     l2 = muColl[0];
     //std::cout<<"test0\n";
     std::map<std::string, bool> ePass = ewp_.CutRecord(l1);
     //std::cout<<"test1\n";
     std::map<std::string, bool> mPass = mwp_.CutRecord(l2,beamSpot);
     //std::cout<<PassAll(mPass)<<std::endl;

     if(scoreHistos_){
      e1Hist_->FillAfterCuts(l1,ePass);
      m1Hist_->FillAfterCuts(l2,beamSpot,mPass);
     }



     if(!(PassAll(ePass)&& PassAll(mPass)))return false;
     TLorentzVector ll1 = Part2LorVec(l1);
     TLorentzVector ll2 = Part2LorVec(l2);
     TLorentzVector l = ll1+ll2;
     if(!((l.M()<=zMassHigh_)&&(l.M()>=zMassLow_)))return false;
     return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
EMuFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EMuFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
EMuFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
EMuFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
EMuFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
EMuFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EMuFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(EMuFilter);
