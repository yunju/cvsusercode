#ifndef patHlttree_h
#define patHlttree_h

#include<iostream>
#include<string>
#include<vector>

#include "TTree.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
class patHltTree{

 public:
  patHltTree(std::string name,TTree* tree,const edm::ParameterSet& iConfig);
  void Fill(const edm::Event& iEvent,const edm::EventSetup& iSetup,const HLTConfigProvider& hltConfig2_);
  void Clear();
 private:
  
  patHltTree(){};
 
  std::string name_;
  TTree* tree_;
  std::vector<int> trigResult_;
  std::vector<int> trigPrescale_;
  std::vector<std::string> trigName_;
  
  std::string HltKeyWordC;

};

#endif
