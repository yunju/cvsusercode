#ifndef __JET_TREE_H_
#define __JET_TREE_H_
/*
  Author: Anil P Singh
   Panjab University
*/





#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DelPanj/TreeMaker/interface/baseTree.h"

using namespace std;
using namespace edm;



class jetTree  : public baseTree{

 public:
  jetTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);//name=patJetAk05
  ~jetTree();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
 
 private:
  jetTree();
  edm::InputTag JetLabel_;
  bool usePFObjects_;
  edm::ParameterSet parSet_;	 

  bool NJet_;
  //Branches common to all the jets.
  std::vector<double> JetPt_;
  std::vector<double> JetEta_;
  std::vector<double> JetPhi_;
  std::vector<double> JetM_;
  std::vector<double> JetRapidity_;
  std::vector<double> JetPx_;
  std::vector<double> JetPy_;
  std::vector<double> JetPz_;
  std::vector<double> JetEn_;
  std::vector<double> JetUnCorrPt_;
  std::vector<double> JetUnCorrPx_;
  std::vector<double> JetUnCorrPy_;
  std::vector<double> JetUnCorrPz_;
  std::vector<double> JetUnCorrEnt_;
  std::vector<double> JetPhotEn_;
  std::vector<double> JetElecEn_;
  std::vector<double> JetMuonEn_;
  std::vector<double> JetHfHadEn_;
  std::vector<double> JetHfEmEn_;
  std::vector<double> JetCharHadE_;
  std::vector<double> JetNeutHadE_;
  std::vector<double> JetCharEmE_;
  std::vector<double> JetCharMuE_;
  std::vector<double> JetNeutEmE_;
  std::vector<double> JetMuonMulti_;
  std::vector<double> JetNeutMulti_;
  std::vector<double> JetCharMulti_;
  std::vector<double> JetCharHadMulti_;
  std::vector<double> JetNeutHadMulti_;
  std::vector<double> JetPhotMulti_;
  std::vector<double> JetElecMulti_;
  std::vector<double> JetHfHadMulti_;
  std::vector<double> JetHfEmMulti_;
  std::vector<double> JetPhotEnFr_;
  std::vector<double> JetMuonEnFr_;
  std::vector<double> JetHfHadEnFr_;
  std::vector<double> JetHfEmEnFr_;
  std::vector<double> JetNeutEmEFr_;
  std::vector<double> JetCharHadEFr_;
  std::vector<double> JetNeutHadEFr_;
  std::vector<double> JetCharEmEFr_;
  std::vector<double> JetCharMuEFr_;
  std::vector<double> JetN60_;
  std::vector<double> JetN90_;
  std::vector<double> JetEmFr_;
  std::vector<double> JetHadFr_;
  std::vector<double> JetEmEbEn_;
  std::vector<double> JetEmEeEn_;
  std::vector<double> JetEmHfEn_;
  std::vector<double> JetHadHbEn_;
  std::vector<double> JetHadHeEn_;
  std::vector<double> JetHadHfEn_;
  std::vector<double> JetHadHoEn_;
  

  std::vector<double> JetTotalEm_; 
  std::vector<double> JetTotalHad_ ;
  std::vector<double> JetEmEbFr_;
  std::vector<double> JetEmEeFr_ ;
  std::vector<double> JetEmHfFr_ ;
  std::vector<double> JetHadHbFr_ ;
  std::vector<double> JetHadHeFr_ ;
  std::vector<double> JetHadHfFr_ ;
  std::vector<double> JetHadHoFr_ ;
  std::vector<double> JetNConstituents_ ; 


};

#endif
