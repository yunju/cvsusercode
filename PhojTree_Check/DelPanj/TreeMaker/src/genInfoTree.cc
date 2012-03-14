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
genInfoTree::AddBranch(int* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x);


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

  edm::Handle<edm::View<pat::Jet> > JetsHandle; 
 
  if( not iEvent.getByLabel("selectedPatJetsPFlow",JetsHandle)){  ///selectedPatJetsAK5PF",JetsHandle)){ 
     edm::LogInfo("GenAnalyzer") << "genJets not found, "
       "skipping event";
         std::cout<<"FATAL EXCEPTION(photon tree  jet info (No selectedPatJetsPFlow)): "<<std::endl; 
     return;
   }

//   edm::Handle<reco::GenJetCollection> genJetsHandle;
  
  edm::Handle<std::vector<reco::GenParticle> >  genJetsHandle;
  if( not iEvent.getByLabel("genParticles",genJetsHandle)){ 
     edm::LogInfo("GenAnalyzer") << "genJets not found, "
       "skipping event";
         std::cout<<"FATAL EXCEPTION(photon tree gen jet info (No GenParticle)): "<<std::endl; 
     return;
   }

   if (JetsHandle.isValid()&&genJetsHandle.isValid())
   {
      for (View<pat::Jet>::const_iterator gjet = JetsHandle->begin(); gjet != JetsHandle->end(); ++gjet) {
      if (!(*gjet).genJet()) continue;
      if ((*gjet).genJet()->pt() < 20) continue;
     
 
      genJetPx_.push_back((*gjet).genJet()->px());
      genJetPy_.push_back((*gjet).genJet()->py()); 
      genJetPz_.push_back((*gjet).genJet()->pz()); 
      genJetE_.push_back((*gjet).genJet()->energy()); 
      genJetP_.push_back((*gjet).genJet()->pt()); 
      genJetPt_.push_back((*gjet).genJet()->pt());
      genJetEta_.push_back((*gjet).genJet()->eta());
      genJetPhi_.push_back((*gjet).genJet()->phi());
      genJetTheta_.push_back((*gjet).genJet()->theta());
      genJetEt_.push_back((*gjet).genJet()->et());

       

     }

   }
 
   ////Get PileupInfo
   edm::InputTag PileupSrc_("addPileupInfo");
   edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
   iEvent.getByLabel(PileupSrc_, PupInfo);
   std::vector<PileupSummaryInfo>::const_iterator PVI;
   for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
     numOfPUInteractions_ = PVI->getPU_NumInteractions();
     //trueNumOfInteractions_ = PVI->getTrueNumInteractions(); ////NOT YET AVAILABLE, BUT SHOULD BE USED TO RE-WEIGHT AS OF CMSSW_4_2_8 ON                                                                                              
   }//end of PVI

}



void  
genInfoTree::SetBranches(){
/*
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
 */ 
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
  AddBranch(&numOfPUInteractions_,"numOfPUInteractions_");
}


void  
genInfoTree::Clear(){
  numOfPUInteractions_ = -99999;

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



