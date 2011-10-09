/*
Lovedeep Kaur Saini
Panjab University, 
Chandigarh
*/
#include <iostream>
#include <cstdlib>
#include "DelPanj/TreeMaker/interface/eventInfo.h"


eventInfo::eventInfo(std::string name, TTree* tree, const edm::ParameterSet& iConfig){
  name_=name;
  tree_=tree;
  
  VertexProducerC       = iConfig.getParameter<edm::InputTag>("VertexProducerPY");
  BeamSpotProducerC     = iConfig.getParameter<edm::InputTag>("BeamSpotProducerPY");
 

  SetBranches();
}


eventInfo::~eventInfo(){
  delete tree_;
}


void
eventInfo::Fill(const edm::Event& iEvent){
  nEvt_   = iEvent.id().event();
  nRun_   = iEvent.id().run();
  nLumiS_ = iEvent.luminosityBlock();
  bunchX_ = iEvent.bunchCrossing();

            // Get the Beam Spot
   reco::BeamSpot beamSpot;
   edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
   iEvent.getByLabel(BeamSpotProducerC,recoBeamSpotHandle);
   beamSpot = *recoBeamSpotHandle;

   // beam spot position
   bspotPos_[0] = beamSpot.x0();
   bspotPos_[1] = beamSpot.y0();
   bspotPos_[2] = beamSpot.z0();


   // Get the primary event vertex
   edm::Handle<reco::VertexCollection> vertexHandle;
   iEvent.getByLabel(VertexProducerC, vertexHandle);

   int nVtxNotFake = 0;

   //  cout <<" start storeVertex "<<endl;
   findGoodVtx          = false;
   
   reco::VertexCollection vertexCollection = *(vertexHandle.product());
   vtx_.SetXYZ(0.,0.,0.);
   double chi2(-1), ndof(-1), normChi2(-1), vtxXError(-1),  vtxYError(-1), vtxZError(-1);
   int vtxNTrk(0), vtxNTrkWeight05(0), nVtxGood(0);
  
   //std::cout<<"check vertexhan"<<vertexHandle.isValid()<<std::endl; 
   for (size_t i=0; i<vertexHandle->size(); ++i)
   {

     //std::cout<<"check vertex why"<<(*vertexHandle)[i].x()<<std::endl;    


   }
   if(vertexCollection.size()>0){
   //std::cout<<"check vertex"<<std::endl; 
 
    reco::VertexCollection::const_iterator vert = vertexCollection.begin();

     for( ; vert!=vertexCollection.end(); ++vert)
     {
       //std::cout<<"check vertex2 "<<vert->x()<<" "<<vert->y()<<" "<<vert->z()<<" "<<vert->ndof()<<std::endl;       
        Double_t vertex_rho = sqrt(vert->x()*vert->x()+
                                  vert->y()*vert->y());

       if ( !(vert->isFake()) && vert->ndof() > 4 && fabs(vert->z()) < 24.0   && vertex_rho < 2.0 )
       {
           //std::cout<<"check vertex3"<<std::endl; 
           if(nVtxGood==0)vtx_ = vert->position();  // set vtx_ to be the first and best vertex
           nVtxGood++;
           findGoodVtx=true;
       }//(!(vert->isFake()) && vert->ndof() > 4 && fabs(vert->z()) < 24.0   && vertex_rho < 2.0)
      

       if ( !(vert->isFake()) )
       {
         //std::cout<<"check vertex4"<<std::endl; 
         vtxXError = vert->xError();
         vtxYError = vert->yError();
         vtxZError = vert->zError();
         chi2      = vert->chi2();
         ndof      = vert->ndof();
         normChi2  = vert->normalizedChi2();
         vtxNTrk   = vert->tracksSize();

         vtxNTrkWeight05 = 0;
         reco::Vertex::trackRef_iterator ittrk;
         for(ittrk = vert->tracks_begin(); ittrk!= vert->tracks_end(); ++ittrk)
           if ( vert->trackWeight(*ittrk) > 0.5 ) vtxNTrkWeight05++;
           
          std::cout<<"check vertex"<<vert->position().X()<<std::endl; 
          vertexX_.push_back(vert->position().X());
          vertexY_.push_back(vert->position().Y());
          vertexZ_.push_back(vert->position().Z());
          vertexXError_.push_back(vtxXError);
          vertexYError_.push_back(vtxYError);
          vertexZError_.push_back(vtxZError);
          vertexChi2_.push_back(chi2);
          vertexNormChi2_.push_back(normChi2);
          vertexNdof_.push_back(ndof);
          nVtxNotFake++;

      }// if ( !(vert->isFake()) )
     } // for( ; vert!=vertexCollection.end(); ++vert)
   }//if (vertexCollection.size()>0)
    nVtxNotFake_=nVtxNotFake;
    nVtxGood_=nVtxGood;
  
 //std::cin.get();
}


void
eventInfo::SetBranches(){
  AddBranch(&nEvt_,"EventNum");
  AddBranch(&nRun_,  "RunNum");
  AddBranch(&nLumiS_, "LumiSection");
  AddBranch(&bunchX_, "BunchXing");
  
  AddBranch(&nVtxGood_, "nVtxGood");
  AddBranch(&nVtxNotFake_, "nVtxNotFake");
  AddBranch(&vertexX_,"vertexX");
  AddBranch(&vertexY_,"vertexY");
  AddBranch(&vertexZ_,"vertexZ");
  AddBranch(&vertexXError_,"vertexXError");
  AddBranch(&vertexYError_,"vertexYError");
  AddBranch(&vertexZError_,"vertexZError");
  AddBranch(&vertexChi2_,"vertexChi2");
  AddBranch(&vertexNormChi2_,"vertexNormChi2");
  AddBranch(&vertexNdof_,"vertexNdof");
  

}

void 
eventInfo::AddBranch(int* x, std::string name){
  std::string brName="EvtInfo_"+name;
  tree_->Branch(brName.c_str(),x,(brName+"/I").c_str());
}
void eventInfo::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName ="EvtInfo_"+name;
  tree_->Branch(brName.c_str(),vec);
}




void 
eventInfo::Clear(){
  nEvt_   = -99999;
  nRun_   = -99999;
  nLumiS_ = -99999;
  bunchX_ = -99999;
  nVtxGood_=-99999;
  nVtxNotFake_ = -99999;

   vertexX_.clear();
   vertexY_.clear();
   vertexZ_.clear();
   vertexXError_.clear();
   vertexYError_.clear();
   vertexZError_.clear();
   vertexChi2_.clear();
   vertexNormChi2_.clear();
   vertexNdof_.clear();




}

