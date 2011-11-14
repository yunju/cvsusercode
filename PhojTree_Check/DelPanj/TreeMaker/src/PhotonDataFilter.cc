// -*- C++ -*-
//
// Original Author:  Anil Pratap Singh,32 2-C17,+41227676591,
//         Created:  Sat Sep 10 18:46:42 CEST 2011
// $Id: PhotonDataFilter.cc,v 1.1 2011/10/09 23:08:46 yunju Exp $
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
   class PhotonDataFilter : public edm::EDFilter {
   public:
     
     explicit PhotonDataFilter(const edm::ParameterSet&);
     ~PhotonDataFilter();
     
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
      double phoTrkIsoMaxC;
      double phoEBSieieMaxC;
      double phoEESieieMaxC;
      



      int Npcrosess;
      int Npass;
    

 
   };
   
   
   PhotonDataFilter::PhotonDataFilter(const edm::ParameterSet& iConfig)
   {
        photonLabel_          = iConfig.getParameter<edm::InputTag> ("photonLabel");
        JetLabel_             = iConfig.getParameter<edm::InputTag> ("jetLabel"); 
        JetetaMaxC            = iConfig.getUntrackedParameter<double>("JetetaMax",5 );
        ptJetMinC             = iConfig.getUntrackedParameter<double>("JetPtMin", 20);
        ptMinC                = iConfig.getUntrackedParameter<double>("GammaPtMin", 40);
        etaMaxC               = iConfig.getUntrackedParameter<double>("GammaEtaMax",2.5);
        hadEmMaxC             = iConfig.getUntrackedParameter<double>("GammaHadEmMax",0.05);
        phoTrkIsoMaxC         = iConfig.getUntrackedParameter<double>("GammaTisoMax",20);
        phoEBSieieMaxC        = iConfig.getUntrackedParameter<double>("GammaEBSieieMax",0.024);
        phoEESieieMaxC        = iConfig.getUntrackedParameter<double>("GammaEESieieMax",0.1); 
        Npcrosess =0;
        Npass =0;
   } 
   
   
   PhotonDataFilter::~PhotonDataFilter()
   {
     
     // do anything here that needs to be done at desctruction time
     // (e.g. close files, deallocate resources etc.)
     
   }
   
   
   //
   // member functions
   //
   
   // ------------ method called on each new Event  ------------
   bool
   PhotonDataFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
          if(fabs(photon.p4().eta())<1.56 && fabs(photon.p4().eta())>1.45  ) continue;
          if(photon.hadronicOverEm()< hadEmMaxC) continue;
          if(fabs(photon.p4().eta()) > etaMaxC) continue;
          if(photon.trkSumPtHollowConeDR04() > phoTrkIsoMaxC) continue;
          if(fabs(photon.p4().eta())<1.45&&photon.sigmaIetaIeta() > phoEBSieieMaxC ) continue;         
          if(fabs(photon.p4().eta())<2.5 && fabs(photon.p4().eta())>1.56 && photon.sigmaIetaIeta() > phoEESieieMaxC ) continue;
          
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
PhotonDataFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonDataFilter::endJob() {
std::cout<<"Npass "<< Npass<<" Npcrosess "<<Npcrosess<<std::endl;
}

// ------------ method called when starting to processes a run  ------------
bool 
PhotonDataFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
PhotonDataFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
PhotonDataFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
PhotonDataFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

int 
PhotonDataFilter::Nev() {
return Npcrosess;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonDataFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PhotonDataFilter);

//  LocalWords:  Brl
