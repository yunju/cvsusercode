// -*- C++ -*-
//
// Original Author:  Anil Pratap Singh,32 2-C17,+41227676591,
//         Created:  Sat Sep 10 18:46:42 CEST 2011
// $Id: PhotonFilter.cc,v 1.4 2012/02/04 04:24:49 yunju Exp $
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

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

//
// class declaration
//
   class PhotonFilter : public edm::EDFilter {
   public:
     
     explicit PhotonFilter(const edm::ParameterSet&);
     ~PhotonFilter();
     
     static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
     
   private:
     virtual void beginJob() ;
     virtual bool filter(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
     
     virtual bool beginRun(edm::Run&, edm::EventSetup const&);
     virtual bool endRun(edm::Run&, edm::EventSetup const&);
     virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
     virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
     virtual int  Nev();
     edm::InputTag photonLabel_;
     edm::InputTag JetLabel_; 
     // ----------member data ---------------------------
      double ptJetMinC;
      double JetetaMaxC;
      double ptMinC;
      double etaMaxC;
      double hadEmMaxC;
      int Npcrosess;
      int Npass;
    

 
   };
   
   
   PhotonFilter::PhotonFilter(const edm::ParameterSet& iConfig)
   {
        photonLabel_          = iConfig.getParameter<edm::InputTag> ("photonLabel");
        JetLabel_             = iConfig.getParameter<edm::InputTag> ("jetLabel"); 
        JetetaMaxC            = iConfig.getUntrackedParameter<double>("JetetaMax",5 );
        ptJetMinC             = iConfig.getUntrackedParameter<double>("JetPtMin", 30);
        ptMinC                = iConfig.getUntrackedParameter<double>("GammaPtMin", 85);
        etaMaxC               = iConfig.getUntrackedParameter<double>("GammaEtaMax",2.6);
        hadEmMaxC             = iConfig.getUntrackedParameter<double>("GammaHadEmMax",0.05);
        Npcrosess =0;
        Npass =0;
   } 
   
   
   PhotonFilter::~PhotonFilter()
   {
     
     // do anything here that needs to be done at desctruction time
     // (e.g. close files, deallocate resources etc.)
     
   }
   
   
   //
   // member functions
   //
   
   // ------------ method called on each new Event  ------------
   bool
   PhotonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
   {
      
      //for checking filter efficiency
      Npcrosess++;
 
    //fetch the input collection
     edm::Handle<pat::PhotonCollection> photonHandle;
     if(not iEvent.getByLabel(photonLabel_,photonHandle)){
      std::cout<<"FATAL EXCEPTION(Photon filter): "<<"Following Not Found: "<<photonLabel_<<std::endl; 
      exit(0);
      }  
       pat::PhotonCollection phColl(*(photonHandle.product()));

       //sort the objects by transverse momentum
 
       
       pat::PhotonCollection::const_iterator ph;
       

       bool foundpho=false;
       for(ph=phColl.begin(); ph!=phColl.end(); ph++)
       {
          pat::Photon photon = pat::Photon(*ph);
          if(photon.pt()< ptMinC ) continue;
          if(photon.hadronicOverEm()> hadEmMaxC) continue;
          if(fabs(photon.p4().eta()) > etaMaxC) continue;
          foundpho = true;
       }
      
      edm::Handle<std::vector<pat::Jet> > JetHandle;
      if(not iEvent.getByLabel(JetLabel_,JetHandle)){
      std::cout<<"FATAL EXCEPTION(Photon filter): "<<"Following Not Found: "
         <<JetLabel_<<std::endl; exit(0);}

        bool foundjet=false;

        const std::vector<pat::Jet>* jets = JetHandle.product();
        std::vector<pat::Jet>::const_iterator jet =jets->begin();
        for(;jet!=jets->end();jet++)
        {
           if(jet->pt()< ptJetMinC ) continue;
           if(fabs(jet->eta())>JetetaMaxC) continue;
           foundjet=true;            
        } 


 
       if(foundpho==true&&foundjet==true){Npass++;}
       //std::cout<<"Npass "<< Npass<<" Npcrosess "<<Npcrosess<<std::endl;
  
       if(!foundpho||!foundjet) return false;
       return true;
   }

// ------------ method called once each job just before starting event loop  ------------
void 
PhotonFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonFilter::endJob() {
std::cout<<"Npass "<< Npass<<" Npcrosess "<<Npcrosess<<std::endl;
}

// ------------ method called when starting to processes a run  ------------
bool 
PhotonFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
PhotonFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
PhotonFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
PhotonFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

int 
PhotonFilter::Nev() {
return Npcrosess;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PhotonFilter);

//  LocalWords:  Brl
