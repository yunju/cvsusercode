#ifndef __eventInfo__
#define __eventInfo__

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h" 
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class eventInfo{

 public:
  eventInfo(std::string name, TTree* tree,const edm::ParameterSet& iConfig);
  ~eventInfo();
  void Fill(const edm::Event& iEvent); 
  void Clear();
 private:
  eventInfo(){};//Don't allow user
  void SetBranches();
  void AddBranch(int* x, std::string name);
  void AddBranch(std::vector<double>* vec, std::string name); 
  edm::InputTag BeamSpotProducerC;
  edm::InputTag VertexProducerC;
  
  TTree *tree_;
  std::string name_;
  int nEvt_;
  int nRun_;
  int nLumiS_;
  int bunchX_;
    
  //vertax information 
  float bspotPos_[3];

  math::XYZPoint vtx_;
  bool    findGoodVtx; 
  int nVtxGood_;
  int nVtxNotFake_;
  std::vector<double> vertexX_;
  std::vector<double> vertexY_;
  std::vector<double> vertexZ_;
  std::vector<double> vertexXError_;
  std::vector<double> vertexYError_;
  std::vector<double> vertexZError_;
  std::vector<double> vertexChi2_;
  std::vector<double> vertexNormChi2_;
  std::vector<double> vertexNdof_;












};

#endif

