#include "DelPanj/TreeMaker/interface/patHltTree.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h" 

patHltTree::patHltTree(std::string name,TTree* tree,const edm::ParameterSet& iConfig)
{
    
  name_=name;
  tree_=tree;

  tree_->Branch("trigResults",&trigResult_);
  tree_->Branch("trigPrescale",&trigPrescale_);
  tree_->Branch("trigName",&trigName_);
 
  HltKeyWordC   = iConfig.getParameter<std::string>("HltKeyWordPY"); 
}

void
patHltTree::Fill(const edm::Event& iEvent,const edm::EventSetup& iSetup,const HLTConfigProvider& hltConfig2_)
{
  Clear();
  using namespace edm;
  
  edm::Handle<edm::TriggerResults> trigResults;

  edm::InputTag trigTag("TriggerResults::HLT");
  if (not iEvent.getByLabel(trigTag, trigResults)) {
    std::cout << ">>> TRIGGER collection does not exist !!!\n";
    return;
  }

  const edm::TriggerNames & trigNames = iEvent.triggerNames(*trigResults);

  for (unsigned int i=0; i<trigResults->size(); i++)
  {
    std::string trigName = trigNames.triggerName(i);
    //reqire to store photon trigger only    
    size_t found;
    std::string hasstr = HltKeyWordC;
    found=trigName.find(hasstr);
    if(found!=std::string::npos)
    {
      //trigName_.push_back(trigName);
      int trigResult = trigResults->accept(i); //bool not to use
     // trigResult_.push_back(trigResult);
      std::string * trigNameTempl = &trigName;  
      //std::cout<<" *trigNameTempl "<<*trigNameTempl<<std::endl;
      //std::cin.get();
      int hltPrescale = hltConfig2_.prescaleValue(iEvent, iSetup, *trigNameTempl);
      bool changed(true);
      //std::cout<<" hlt "<<trigName <<" "<<hltPrescale<<"  "<<trigResult <<std::endl;
    //std::cin.get();
      if (hltPrescale!=1||trigResult!=1) continue;
      trigName_.push_back(trigName);      
      trigResult_.push_back(trigResult);
      trigPrescale_.push_back(hltPrescale);
      //std::cout<<" hltpre1 "<<trigName <<" "<<hltPrescale<<"  "<<trigResult<<std::endl;
      //std::cin.get();    



    } 
 
  }

    /* 
 for ( std::vector<std::string>::const_iterator trigNameTempl = ttrigResults.begin();  trigNameTempl != ttrigResults.end(); ++ttrigResults) {

    std::cout<<" strigNameTempl "<< hltConfig_.prescaleValue(iEvent, iSetup, *trigNameTempl)<<std::endl;
  }

*/

}


void
patHltTree::Clear(){
  trigResult_.clear();
  trigName_.clear();
}


