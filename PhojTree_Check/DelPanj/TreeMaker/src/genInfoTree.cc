/*
Anil Singh
Panjab University

*/

#include "DelPanj/TreeMaker/interface/genInfoTree.h"

//---------------------------------------------------------------
// Add Branches to the genTree
//---------------------------------------------------------------
void
genInfoTree::AddBranch(double* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,(brName+"/D").c_str());

}


//---------------------------------------------------------------
//---------------------------------------------------------------
void
genInfoTree::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);


}

//---------------------------------------------------------------
//---------------------------------------------------------------
void
genInfoTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

genInfoTree::genInfoTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig)
{
  tree_=tree; 
  genPartLabel_ = iConfig.getParameter<edm::InputTag>("genPartLabel");
  SetBranches();
}


genInfoTree::~genInfoTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete tree_;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
genInfoTree::Fill(const edm::Event& iEvent)
{
  Clear();
   using namespace edm;/*
   edm::Handle<reco::GenParticleCollection> genParticleHandle;
   if(not iEvent.getByLabel(genPartLabel_, genParticleHandle))
     {
       std::cout<<
	 "GenAnalyzer: Generator Level Information not found\n"
		<<std::endl;
     }


   const reco::GenParticleCollection* genColl= &(*genParticleHandle);
   reco::GenParticleCollection::const_iterator geni = genColl->begin();
   for(; geni!=genColl->end();geni++){
     reco::GenParticle gen = *geni;
     
     //Look out for the GenMuons
     if(gen.status()==1){
    
       genParPx_.push_back(gen.px());
       genParPy_.push_back(gen.py());
       genParPz_.push_back(gen.pz());
       genParE_.push_back(gen.energy());
       genParP_.push_back(gen.p());
       genParPt_.push_back(gen.pt());
       genParEta_.push_back(gen.eta());
       genParPhi_.push_back(gen.phi());
       genParTheta_.push_back(gen.theta());
       genParEt_.push_back(gen.et());
       genParQ_.push_back(gen.charge());
       genParId_.push_back(gen.pdgId());
       genParSt_.push_back(gen.status());
     
     }
   }
*/
   edm::Handle<reco::GenJetCollection> genJetsHandle;
   if( not iEvent.getByLabel("iterativeCone5GenJets",genJetsHandle)){ 
     edm::LogInfo("GenAnalyzer") << "genJets not found, "
       "skipping event"; 
     return;
   }
   const reco::GenJetCollection* genJetColl = &(*genJetsHandle);
   reco::GenJetCollection::const_iterator gjeti = genJetColl->begin();
   
   for(; gjeti!=genJetColl->end();gjeti++){
     reco::GenParticle gjet = *gjeti;
     genJetPx_.push_back(gjet.px());
     genJetPy_.push_back(gjet.py()); 
     genJetPz_.push_back(gjet.pz()); 
     genJetE_.push_back(gjet.energy()); 
     genJetP_.push_back(gjet.p()); 
     genJetPt_.push_back(gjet.pt());
     genJetEta_.push_back(gjet.eta());
     genJetPhi_.push_back(gjet.phi());
     genJetTheta_.push_back(gjet.theta());
     genJetEt_.push_back(gjet.et());
   }
   
}



void  
genInfoTree::SetBranches(){
  AddBranch(&genParPx_, "genParPx_");
  AddBranch(&genParPy_,"genParPy_");
  AddBranch(&genParPz_,"genParPz_");
  AddBranch(&genParE_, "genParE_");
  AddBranch(&genParP_,"genParP_");
  AddBranch(&genParTheta_,"genParTheta_");
  AddBranch(&genParPt_, "genParPt_");
  AddBranch(&genParEta_,"genParEta_");
  AddBranch(&genParPhi_,"genParPhi_");
  AddBranch(&genParEt_,"genParEt_");
  AddBranch(&genParQ_,"genParQ_");
  AddBranch(&genParId_,"genParId_");
  AddBranch(&genParSt_,"genParSt_");
  AddBranch(&genJetPx_, "genJetPx_");
  AddBranch(&genJetPy_,"genJetPy_");
  AddBranch(&genJetPz_,"genJetPz_");
  AddBranch(&genJetE_, "genJetE_");
  AddBranch(&genJetP_,"genJetP_");
  AddBranch(&genJetTheta_,"genJetTheta_");
  AddBranch(&genJetPt_,"genJetPt_");
  AddBranch(&genJetEta_,"genJetEta_");
  AddBranch(&genJetPhi_,"genJetPhi_");
  AddBranch(&genJetEt_,"genJetEt_"); 
  AddBranch(&genJetQ_,"genJetQ_");

}


void  
genInfoTree::Clear(){
   genParPx_.clear();
   genParPy_.clear();
   genParPz_.clear();
   genParE_.clear();
   genParP_.clear();
   genParTheta_.clear();
   genParPt_.clear();
   genParEta_.clear();
   genParPhi_.clear();
   genParEt_.clear();
   genParQ_.clear();
   genParId_.clear();
   genParSt_.clear();
   genJetPx_.clear();
   genJetPy_.clear();
   genJetPz_.clear();
   genJetE_.clear();
   genJetP_.clear();
   genJetTheta_.clear();
   genJetPt_.clear();
   genJetEta_.clear(); 
   genJetPhi_.clear(); 
   genJetEt_.clear();
   genJetQ_.clear();
}



