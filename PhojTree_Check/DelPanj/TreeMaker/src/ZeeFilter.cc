// -*- C++ -*-
//
// Original Author:  Anil Pratap Singh,32 2-C17,+41227676591,
//         Created:  Sat Sep 10 18:46:42 CEST 2011
// $Id: ZeeFilter.cc,v 1.4 2011/09/17 01:16:47 lovedeep Exp $
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
#include "DataFormats/Math/interface/LorentzVector.h"

//
// class declaration
//
/*------*/
 bool PassAll(std::map<std::string, bool> cutrecd){
   std::map<std::string, bool>::iterator iter= cutrecd.begin();
   bool decision =1 ;
   for(;iter!=cutrecd.end();iter++){
     decision = decision&&iter->second;     
   }
   return decision;
 }


 class eWkPt{
   /*
    I will consider moving this electron selector
    into an independent file, so that other 
    analyzers may be able use it. Let us use 
    and evolve it here atfirst....  
   */


 public:
   eWkPt(const edm::ParameterSet ps):
     //kinematic parameters
     ptx_  (ps.getParameter<double>("ptx")),
     etax_ (ps.getParameter<double>("etax")),
     //identification parameters
     idBrl_(ps.getParameter<edm::ParameterSet> ("idBrl")),   
     idEcp_(ps.getParameter<edm::ParameterSet> ("idEcp")),
     //isolation parameters 
     usePfIso_(ps.getParameter<bool>("usePfIso")),  
     useRelIso_(ps.getParameter<bool>("useRelIso")),
     useCombIso_(ps.getParameter<bool>("useCombIso")),
     coneRadius_(ps.getParameter<double>("coneRad")),
     isoBrl_(ps.getParameter<edm::ParameterSet> ("isoBrl")),
     isoEcp_(ps.getParameter<edm::ParameterSet> ("isoEcp")){     
     
     
     sieieBrlX_= (idBrl_.getParameter<double>("sieie"));
     delphiBrlX_=(idBrl_.getParameter<double>("delphi"));
     detainBrlX_=(idBrl_.getParameter<double>("detain"));
     distBrlX_=(idBrl_.getParameter<double>("dist"));
     dcotBrlX_=(idBrl_.getParameter<double>("dcot"));
     hoeBrlX_=(idBrl_.getParameter<double>("hoe"));
     nmisHitBrlX_=(idBrl_.getParameter<double>("nmisHit"));
     
     sieieEcpX_= (idEcp_.getParameter<double>("sieie"));
     delphiEcpX_=(idEcp_.getParameter<double>("delphi"));
     detainEcpX_=(idEcp_.getParameter<double>("detain"));
     distEcpX_=(idEcp_.getParameter<double>("dist"));
     dcotEcpX_=(idEcp_.getParameter<double>("dcot"));
     hoeEcpX_=(idEcp_.getParameter<double>("hoe"));
     nmisHitEcpX_=(idEcp_.getParameter<double>("nmisHit"));
   
     if(!useCombIso_){
     iso1BrlX_ =(isoBrl_.getParameter<double>("iso1"));
     iso2BrlX_ =(isoBrl_.getParameter<double>("iso2"));
     iso3BrlX_ =(isoBrl_.getParameter<double>("iso3"));
     }else 
     iso4BrlX_ =(isoBrl_.getParameter<double>("iso4"));
     if(!useCombIso_){
     iso1EcpX_ =(isoEcp_.getParameter<double>("iso1"));
     iso2EcpX_ =(isoEcp_.getParameter<double>("iso2"));
     iso3EcpX_ =(isoEcp_.getParameter<double>("iso3"));
     }else
     iso4EcpX_ =(isoEcp_.getParameter<double>("iso4"));
     
   }
   
   std::map<std::string, bool> CutRecord(pat::Electron& e){
       std::map<std::string, bool> cuts;
       
       //good idea to call the eta and pt methods for once.
       double pt = e.pt();
       double eta = e.eta();
       
       /*

       % ONE OF THE TWO MUST BE TRUE: pflowIso or caloIso
       % ONE OF THESE MUST BE TRUE: Rel         
         --IF USE PFLOW ISO = TURE.
          iso1  = chHadIso
          iso2  = neHadIso
          iso3  = gammaIso
          iso4  = combIso
           Note: I dont know what they are when pat::Electron comes
                 from a standard gsfElelctron. So we will place a
                 check for candiadate being pflow, and will allow
                 these only it has a valid pfcandref embedded.    

        --IF USE CALO ISO = TRUE
         iso1  = trkIso
         iso2  = ecalIso
         iso3  = hcalIso
         iso4  = combIso
           Note: These can be calculated for both PF or Calo
                 electrons.
       */

      //Default values must fail any sensible IsoCut
      double iso1 = 999.;
      double iso2 = 999.;
      double iso3 = 999.;
      double iso4 = 999.;

      if(usePfIso_){   
       iso1  = e.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin( coneRadius_ ).first;
       iso2  = e.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin( coneRadius_ ).first;
       iso3  = e.isoDeposit(pat::PfGammaIso)->depositAndCountWithin( coneRadius_ ).first;
      }
      else {
	iso1   = e.trackIso();
	iso2   = e.ecalIso(); 
	iso3   = e.hcalIso(); 
      }

      iso4  = (iso1+iso2+iso3);
      if(useRelIso_){
	iso1 = iso1/pt;
	iso2 = iso2/pt;
	iso3 = iso3/pt;
	iso4 = iso4/pt;
      }
 
       if(e.isEB()){
	 cuts["ptx"]    = pt>ptx_;
	 cuts["etax"]   = fabs(eta)<etax_;
	 cuts["sieie"]  = fabs(e.scSigmaIEtaIEta())<sieieBrlX_;
	 cuts["delphi"] = fabs(e.deltaPhiSuperClusterTrackAtVtx())<delphiBrlX_;
	 cuts["detain"] = fabs(e.deltaEtaSuperClusterTrackAtVtx())<detainBrlX_;

 	 cuts["dcot"]   = fabs(e.convDcot())> dcotBrlX_;
         cuts["dist"]   = fabs(e.convDist())> distBrlX_;	
         cuts["hoe"]    = e.hadronicOverEm() <hoeBrlX_;
         cuts["nmshits"]= e.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits()<=nmisHitBrlX_;

	 if(useCombIso_){
           cuts["iso1"]   = 1;
           cuts["iso2"]   = 1;
           cuts["iso3"]   = 1;
           cuts["iso4"]   = iso4 < iso4EcpX_;
         }
         else{
           cuts["iso1"]   = iso1 < iso1EcpX_;
           cuts["iso2"]   = iso2 < iso2EcpX_;
           cuts["iso3"]   = iso3 < iso3EcpX_;
           cuts["iso4"]   = 1;
         }

       }
       
       else if(e.isEE()){
	 cuts["ptx"]    = pt>ptx_;
	 cuts["etax"]   = fabs(e.eta())<etax_;
	 cuts["sieie"]  = fabs(e.scSigmaIEtaIEta())<sieieEcpX_;
	 cuts["delphi"] = fabs(e.deltaPhiSuperClusterTrackAtVtx())<delphiEcpX_;
	 cuts["detain"] = fabs(e.deltaEtaSuperClusterTrackAtVtx())<detainEcpX_;
         cuts["dcot"]   = fabs(e.convDcot()) > dcotEcpX_;
         cuts["dist"]   = fabs(e.convDist()) > distEcpX_;
         cuts["hoe"]    = e.hadronicOverEm() < hoeEcpX_;
         cuts["nmshits"]= e.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits()<=nmisHitEcpX_;
	 if(useCombIso_){
	   cuts["iso1"]   = 1;
	   cuts["iso2"]   = 1;
	   cuts["iso3"]   = 1;
	   cuts["iso4"]   = iso4 < iso4EcpX_;
	 }
	 else{
	   cuts["iso1"]   = iso1 < iso1EcpX_;
	   cuts["iso2"]   = iso2 < iso2EcpX_;
	   cuts["iso3"]   = iso3 < iso3EcpX_;
	   cuts["iso4"]   = 1;
	 }
       }
       
       else{
 	 cuts["ptx"]    = 0;
	 cuts["etax"]   = 0;
	 cuts["sieie"]  = 0;
	 cuts["delphi"] = 0;
	 cuts["detain"] = 0;
  	 cuts["dcot"] = 0;
  	 cuts["dist"] = 0;	
         cuts["hoe"]    = 0;
         cuts["nmshits"] = 0;
	 cuts["iso1"]   = 0;
	 cuts["iso2"]   = 0;
	 cuts["iso3"]   = 0;
         cuts["iso4"]   = 0;
        
       }
   
       return cuts;
     }
   


   ~eWkPt(){}
   double ptx_;
   double etax_;
   
   double sieieBrlX_;
   double delphiBrlX_;
   double detainBrlX_;
   double distBrlX_;
   double dcotBrlX_;
   double hoeBrlX_;
   double nmisHitBrlX_;
   double sieieEcpX_;
   double delphiEcpX_;
   double detainEcpX_;
   double distEcpX_;
   double dcotEcpX_;
   double hoeEcpX_;
   double nmisHitEcpX_;

   double iso1BrlX_;
   double iso2BrlX_;
   double iso3BrlX_;
   double iso4BrlX_;

   double iso1EcpX_;
   double iso2EcpX_;
   double iso3EcpX_;
   double iso4EcpX_;



   //some flags.
   bool usePfIso_;
   bool useRelIso_;
   bool useCombIso_;
   bool coneRadius_;

  

   //some pset declaration
   edm::ParameterSet idBrl_;
   edm::ParameterSet idEcp_;
   edm::ParameterSet isoBrl_;
   edm::ParameterSet isoEcp_; 
   };
     
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
			
   
     //electron selection cuts.
     eWkPt ewp1_;//for lead ele.
     eWkPt ewp2_;//for sublead ele
   };
   
   
   ZeeFilter::ZeeFilter(const edm::ParameterSet& iConfig):
     patElecLabel_(iConfig.getParameter<edm::InputTag>("patElectrons")),
     ewp1_(iConfig.getParameter<edm::ParameterSet>("leadElecPset_")),
     ewp2_(iConfig.getParameter<edm::ParameterSet>("subLeadElecPset_"))
     {  }
   
   
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
    
     edm::Handle<pat::ElectronCollection> patElecHandle;
     if(not iEvent.getByLabel(patElecLabel_,patElecHandle)){
       std::cout<<"ZEE filter FATAL EXCEPTION: "<<"Following Not Found: "
		<<patElecLabel_<<std::endl; return false;}
     pat::ElectronCollection eleColl(*(patElecHandle.product()));
     if(eleColl.size()<2)return false;
     std::map<std::string, bool> leadPass = ewp1_.CutRecord(eleColl[0]);
     std::map<std::string, bool> subLeadPass = ewp2_.CutRecord(eleColl[1]);
     if(!(PassAll(leadPass)&& PassAll(subLeadPass)))return false;
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
