#ifndef _GEN_INFO_TREE_H_
#define _GEN_INFO_TREE_H_

#include <memory>
#include <vector>
#include <string>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

//
// class declaration
//

class genInfoTree{

 public:
  genInfoTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~genInfoTree();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
  
  
  // ----------member data ---------------------------
  void AddBranch(double* x, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<int>*, std::string name);
  
  edm::InputTag genPartLabel_;
  
  
  TTree* tree_;

  std::vector<double> genParPx_;
  std::vector<double> genParPy_;
  std::vector<double> genParPz_;
  std::vector<double> genParE_;
  std::vector<double> genParP_;
  std::vector<double> genParTheta_;
  std::vector<double> genParPt_;
  std::vector<double> genParEta_;
  std::vector<double> genParPhi_;
  std::vector<double> genParEt_;
  std::vector<double> genParQ_;
  std::vector<double> genParId_;
  std::vector<double> genParSt_;

  std::vector<double> genJetPx_;
  std::vector<double> genJetPy_;
  std::vector<double> genJetPz_;
  std::vector<double> genJetE_;
  std::vector<double> genJetP_;
  std::vector<double> genJetTheta_;
  std::vector<double> genJetPt_;
  std::vector<double> genJetEta_; 
  std::vector<double> genJetPhi_; 
  std::vector<double> genJetEt_;
  std::vector<double> genJetQ_;
};

#endif
