/*
Lovedeep Kaur Saini
Panjab University, 
Chandigarh
*/
#include <iostream>
#include "DelPanj/TreeMaker/interface/eventInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

eventInfo::eventInfo(std::string name, TTree* tree){
  name_=name;
  tree_=tree;
  SetBranches();

}


eventInfo::~eventInfo(){
  delete tree_;
}


void
eventInfo::Fill(const edm::Event& iEvent){
  Clear();

  nEvt_   = iEvent.id().event();
  nRun_   = iEvent.id().run();
  nLumiS_ = iEvent.luminosityBlock();
  bunchX_ = iEvent.bunchCrossing();
  
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("offlinePrimaryVertices",recVtxs);
  for(unsigned int ind=0;ind<recVtxs->size();ind++) {
    if (!((*recVtxs)[ind].isFake()) && ((*recVtxs)[ind].ndof()>4) 
	&& (fabs((*recVtxs)[ind].z())<=24.0) &&  
	((*recVtxs)[ind].position().Rho()<=2.0) ) {   
      nVtx_++;
      vertexX_.push_back((*recVtxs)[ind].x());
      vertexY_.push_back((*recVtxs)[ind].y());
      vertexZ_.push_back((*recVtxs)[ind].z());
    }
  }
}

void
eventInfo::SetBranches(){
  AddBranch(&nEvt_,"EventNum");
  AddBranch(&nRun_,  "RunNum");
  AddBranch(&nLumiS_, "LumiSection");
  AddBranch(&bunchX_, "BunchXing");
  AddBranch(&nVtx_, "NumVtx");
  AddBranch(&vertexX_, "VertexX_");
  AddBranch(&vertexY_, "VertexY_");
  AddBranch(&vertexZ_, "VertexZ_");
}

void
eventInfo::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName = "EvtInfo_"+name;
  tree_->Branch(brName.c_str(),vec);
}


void 
eventInfo::AddBranch(int* x, std::string name){
  std::string brName="EvtInfo_"+name;
  tree_->Branch(brName.c_str(),x,(brName+"/I").c_str());
}

void 
eventInfo::Clear(){
  nEvt_   = -99999;
  nRun_   = -99999;
  nLumiS_ = -99999;
  bunchX_ = -99999;
  nVtx_ = 0;
  vertexX_.clear();
  vertexY_.clear();
  vertexZ_.clear();
}

