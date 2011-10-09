// -*- C++ -*-
//
// Original Author:  Anil Pratap Singh,32 2-C17,+41227676591,
//         Created:  Sat Sep 10 18:46:42 CEST 2011
// $Id: ZeeFilter.cc,v 1.1 2011/10/05 22:15:29 anil Exp $
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

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/eSelector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "DelPanj/TreeMaker/interface/eHist.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


   class ZeeFilter : public edm::EDFilter {
   public:
     
     explicit ZeeFilter(const edm::ParameterSet&);
     ~ZeeFilter();
     
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
     bool scoreHistos_;  
     int count	; //debugging help		
   
   
     //electron selection cuts.
     eSelector ewp1_;//for lead ele.
     eSelector ewp2_;//for sublead ele
     double zMassLow_;
     double zMassHigh_;
     eHist* e1Hist_;
     eHist* e2Hist_;
   	
   };
   
   
   ZeeFilter::ZeeFilter(const edm::ParameterSet& iConfig):
     patElecLabel_(iConfig.getParameter<edm::InputTag>("patElectrons")),
     ewp1_(iConfig.getParameter<edm::ParameterSet>("leadElecPset_")),
     ewp2_(iConfig.getParameter<edm::ParameterSet>("subLeadElecPset_")),
     zMassLow_(iConfig.getParameter<double>("zMassLowerLimit")),
     zMassHigh_(iConfig.getParameter<double>("zMassUpperLimit")),
     scoreHistos_(iConfig.getParameter<double>("scoreHistos"))
     {
       count=0; //debugging
      edm::Service<TFileService> fs;
      if(scoreHistos_){
      e1Hist_ = new eHist("LeadingElectron", fs);
      e2Hist_ = new eHist("TrailingElectron", fs);
     }
   }
   
   ZeeFilter::~ZeeFilter()
   {
     
     // do anything here that needs to be done at desctruction time
     // (e.g. close files, deallocate resources etc.)
     
   }
   
   
   //
   // member functions
   //
   
   // ------------ method called on each new Event  ------------
   bool
   ZeeFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
   {
     //     count++;
     //std::cout<<"+++++++++++++++++++++++++++\n"<<count<<"\n";
     edm::Handle<pat::ElectronCollection> patElecHandle;
     if(not iEvent.getByLabel(patElecLabel_,patElecHandle)){
       std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
		<<patElecLabel_<<std::endl; return false;}
     pat::ElectronCollection eleColl(*(patElecHandle.product()));
     if(eleColl.size()<2)return false;
     //count++;
     //std::cout<<count<<std::endl;	
     pat::Electron e1 = eleColl[0];
     pat::Electron e2 = eleColl[1];
     std::map<std::string, bool> leadPass = ewp1_.CutRecord(e1);
     std::map<std::string, bool> subLeadPass = ewp2_.CutRecord(e2);
     if(scoreHistos_){
       e1Hist_->FillAfterCuts(e1,leadPass);
       e2Hist_->FillAfterCuts(e2,subLeadPass);
     }
    
     if(!(PassAll(leadPass)&& PassAll(subLeadPass)))return false;
     count++;
     //std::cout<<count<<std::endl;
     TLorentzVector l1 = Part2LorVec(e1);
     TLorentzVector l2 = Part2LorVec(e2);
     TLorentzVector l = l1+l2;
     if(!((l.M()<=zMassHigh_)&&(l.M()>=zMassLow_)))return false;
     //count++;
     //std::cout<<ewp1_.distBrlX_<<"\t"<<count<<std::endl;
     //std::cout<<ewp1_.etaX_<<"\t"<<count<<std::endl;
   
     return true;
   }

// ------------ method called once each job just before starting event loop  ------------
void 
ZeeFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZeeFilter::endJob() {
//std::cout<<count<<std::endl;	
//e1Hist_->SetProperRanges();
}

// ------------ method called when starting to processes a run  ------------
bool 
ZeeFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
ZeeFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
ZeeFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
ZeeFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZeeFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ZeeFilter);

//  LocalWords:  Brl
