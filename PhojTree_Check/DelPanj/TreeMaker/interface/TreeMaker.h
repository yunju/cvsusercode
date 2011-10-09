#ifndef  TREE_MAKER_H
#define  TREE_MAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DelPanj/TreeMaker/interface/eventInfo.h"
#include "DelPanj/TreeMaker/interface/genInfoTree.h"
#include "DelPanj/TreeMaker/interface/patMuonTree.h"
#include "DelPanj/TreeMaker/interface/patElecTree.h"
#include "DelPanj/TreeMaker/interface/patMetTree.h"
#include "DelPanj/TreeMaker/interface/patHltTree.h"
#include "DelPanj/TreeMaker/interface/jetTree.h"
#include "DelPanj/TreeMaker/interface/photonTree.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


#include "TTree.h"
#include "TFile.h"


class TreeMaker : public edm::EDAnalyzer {
   public:
      explicit TreeMaker(const edm::ParameterSet&);
      ~TreeMaker();


   private:
      virtual void beginJob() ;
      virtual void beginRun(const edm::Run & iRun, const edm::EventSetup & iSetup);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      TFile* file;
      TTree* tree_;
      std::string outFileName_ ;

      HLTConfigProvider hltConfig_; 
      int NprocessTree;      
 
      bool fillEventInfo_;
      bool fillGenInfo_;
      bool fillTrigInfo_;
      bool fillMuonInfo_;
      bool fillElecInfo_;
      bool fillMetInfo_;
      bool fillJetInfo_;
      bool fillPhotInfo_;
      
      
      eventInfo   *eventInfo_;
      genInfoTree *genInfoTree_;    
      patMuonTree *patMuTree_;
      patElecTree *patElecTree_;
      patMetTree  *patMetTree_;
      jetTree     *jetTree_;
      photonTree  *photonTree_;
      patHltTree  *patHltTree_;      

};


#endif
