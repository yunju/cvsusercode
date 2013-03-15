// -*- C++ -*-
//
// Package:    AnalyzerCheck
// Class:      AnalyzerCheck
// 
/**\class AnalyzerCheck AnalyzerCheck.cc AnaCheck/AnalyzerCheck/src/AnalyzerCheck.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yun-Ju Lu,27 2-004,+41227676186,
//         Created:  Fri Mar 15 11:27:56 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
//
// class declaration
//

class AnalyzerCheck : public edm::EDAnalyzer {
   public:
      explicit AnalyzerCheck(const edm::ParameterSet&);
      ~AnalyzerCheck();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
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
AnalyzerCheck::AnalyzerCheck(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


AnalyzerCheck::~AnalyzerCheck()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
AnalyzerCheck::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<std::vector<pat::Electron> > patElecHandle;
   iEvent.getByLabel("userDataSelectedElectrons",patElecHandle); 
   



}


// ------------ method called once each job just before starting event loop  ------------
void 
AnalyzerCheck::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AnalyzerCheck::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
AnalyzerCheck::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
AnalyzerCheck::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
AnalyzerCheck::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
AnalyzerCheck::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AnalyzerCheck::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzerCheck);
