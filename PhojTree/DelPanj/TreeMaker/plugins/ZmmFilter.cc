// -*- C++ -*-
//
// Original Author:  Anil Pratap Singh,32 2-C17,+41227676591,
//         Created:  Sat Sep 10 18:46:42 CEST 2011
// $Id: ZmmFilter.cc,v 1.1 2011/10/05 22:15:29 anil Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <algorithm>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DelPanj/TreeMaker/interface/muHist.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/muSelector.h"
#include "DelPanj/TreeMaker/interface/utils.h"

   class ZmmFilter : public edm::EDFilter {
   public:
     
     explicit ZmmFilter(const edm::ParameterSet&);
     ~ZmmFilter();
     
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
  
     //mutron tag
     edm::InputTag patMuLabel_;
     edm::InputTag beamSpotLabel_;
     bool scoreHistos_;			
   
     //muon smution cuts.
     muSelector mwp1_;//for lead ele.
     muSelector mwp2_;//for sublead ele
     double zMassLow_;
     double zMassHigh_;
     muHist* m1Hist_;
     muHist* m2Hist_;
   };
   
   
   ZmmFilter::ZmmFilter(const edm::ParameterSet& iConfig):
     patMuLabel_(iConfig.getParameter<edm::InputTag>("patMuons")),
     beamSpotLabel_(iConfig.getParameter<edm::InputTag>("beamSpot")),
     mwp1_(iConfig.getParameter<edm::ParameterSet>("leadMuPset_")),
     mwp2_(iConfig.getParameter<edm::ParameterSet>("subLeadMuPset_")),
     zMassLow_(iConfig.getParameter<double>("zMassLowerLimit")),
     zMassHigh_(iConfig.getParameter<double>("zMassUpperLimit")),
     scoreHistos_(iConfig.getParameter<double>("scoreHistos"))
     { 

      edm::Service<TFileService> fs;
      if(scoreHistos_){
      m1Hist_ = new muHist("leadingMu", fs);
      m2Hist_ = new muHist("subLeadingMu", fs);
      }

     }
   
   ZmmFilter::~ZmmFilter()
   {
     
     // do anything here that needs to be done at desctruction time
     // (e.g. close files, deallocate resources etc.)
     
   }
   
   
   //
   // member functions
   //
   
   // ------------ method called on each new Event  ------------
   bool
   ZmmFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
   {

     //add beamspot also.
     edm::Handle<reco::BeamSpot> beamSpotHandle;
     iEvent.getByLabel(beamSpotLabel_, beamSpotHandle);
     reco::BeamSpot beamSpot = (*beamSpotHandle);

     edm::Handle<pat::MuonCollection> patMuHandle;
     if(not iEvent.getByLabel(patMuLabel_,patMuHandle)){
       std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
		<<patMuLabel_<<std::endl; return false;}
     pat::MuonCollection muColl(*(patMuHandle.product()));
     if(muColl.size()<2)return false;
     pat::Muon e1 = muColl[0];
     pat::Muon e2 = muColl[1];

     std::map<std::string, bool> leadPass = mwp1_.CutRecord(e1,beamSpot);
     std::map<std::string, bool> subLeadPass = mwp2_.CutRecord(e2,beamSpot);
     if(scoreHistos_){
       m1Hist_->FillAfterCuts(e1,beamSpot,leadPass);
       m2Hist_->FillAfterCuts(e2,beamSpot,subLeadPass);
     } 
     if(!(PassAll(leadPass)&& PassAll(subLeadPass)))return false;
     TLorentzVector l1 = Part2LorVec(e1);
     TLorentzVector l2 = Part2LorVec(e2);
     TLorentzVector l = l1+l2;
     if(!((l.M()<=zMassHigh_)&&(l.M()>=zMassLow_)))return false;
     return true;
   }

// ------------ method called once each job just before starting event loop  ------------
void 
ZmmFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZmmFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
ZmmFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
ZmmFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
ZmmFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
ZmmFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZmmFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ZmmFilter);

//  LocalWords:  Brl
