#include "DelPanj/TreeMaker/interface/eSelector.h"

eSelector::eSelector(const edm::ParameterSet ps):      
  ptX_  (ps.getParameter<double>("ptx")),
  etaX_ (ps.getParameter<double>("etax")),
  //identification parameters
  idBrl_(ps.getParameter<edm::ParameterSet> ("idBrl")),   
  idEcp_(ps.getParameter<edm::ParameterSet> ("idEcp")),
  //isolation parameters 
  usePfIso_(ps.getParameter<bool>("usePfIso")),  
  useRelIso_(ps.getParameter<bool>("useRelIso")),
  useCombIso_(ps.getParameter<bool>("useCombIso")),
  coneRadius_(ps.getParameter<double>("coneRad")),
  isoBrl_(ps.getParameter<edm::ParameterSet> ("isoBrl")),
  isoEcp_(ps.getParameter<edm::ParameterSet> ("isoEcp"))

{     
  
  sieieBrlX_= (idBrl_.getParameter<double>("sieie"));
  delphiBrlX_=(idBrl_.getParameter<double>("delphi"));
  detainBrlX_=(idBrl_.getParameter<double>("detain"));
  distBrlX_=(idBrl_.getParameter<double>("dist"));
  dcotBrlX_= idBrl_.getParameter<double>("dcot");
  hoeBrlX_=(idBrl_.getParameter<double>("hoe"));
  nmisHitBrlX_=(idBrl_.getParameter<double>("nmisHit"));
 

//std::cout<<"iiiiiiiiiiiiiiiiiioooo"<<idBrl_.getParameter<double>("dist")<<std::endl;


 
  sieieEcpX_= (idEcp_.getParameter<double>("sieie"));
  delphiEcpX_=(idEcp_.getParameter<double>("delphi"));
  detainEcpX_=(idEcp_.getParameter<double>("detain"));
  distEcpX_=(idEcp_.getParameter<double>("dist"));
  dcotEcpX_=(idEcp_.getParameter<double>("dcot"));
  hoeEcpX_=(idEcp_.getParameter<double>("hoe"));
  nmisHitEcpX_=(idEcp_.getParameter<double>("nmisHit"));
  
  if(!useCombIso_){
    iso1BrlX_ =(isoBrl_.getParameter<double>("iso1"));
    iso2BrlX_ =(isoBrl_.getParameter<double>("iso2"));
    iso3BrlX_ =(isoBrl_.getParameter<double>("iso3"));
  }else 
    iso4BrlX_ =(isoBrl_.getParameter<double>("iso4"));
  if(!useCombIso_){
    iso1EcpX_ =(isoEcp_.getParameter<double>("iso1"));
    iso2EcpX_ =(isoEcp_.getParameter<double>("iso2"));
    iso3EcpX_ =(isoEcp_.getParameter<double>("iso3"));
  }else
    iso4EcpX_ =(isoEcp_.getParameter<double>("iso4"));
  
}


std::map<std::string, bool> 
eSelector::CutRecord(pat::Electron& e){
       std::map<std::string, bool> cuts;
       
       //good idea to call the eta and pt methods for once.
       double pt = e.pt();
       double eta = e.eta();
       bool ingap= fabs(eta)>1.446 && fabs(eta)< 1.566;
//       std::cout<<"isoBrlX: "<< dcotBrlX_<<std::endl;
       /*

       % ONE OF THE TWO MUST BE TRUE: pflowIso or caloIso
       % ONE OF THESE MUST BE TRUE: Rel         
         --IF USE PFLOW ISO = TURE.
          iso1  = chHadIso
          iso2  = neHadIso
          iso3  = gammaIso
          iso4  = combIso
           Note: I dont know what they are when pat::Electron comes
                 from a standard gsfElelctron. So we will place a
                 check for candiadate being pflow, and will allow
                 these only it has a valid pfcandref embedded.    

        --IF USE CALO ISO = TRUE
         iso1  = trkIso
         iso2  = ecalIso
         iso3  = hcalIso
         iso4  = combIso
           Note: These can be calculated for both PF or Calo
                 electrons.
       */
      if(!ingap){
      //Default values must fail any sensible IsoCut
      double iso1 = 999.;
      double iso2 = 999.;
      double iso3 = 999.;
      double iso4 = 999.;

      if(usePfIso_){  
  //    std::cout<<"singRelIso"<<std::endl;
	iso1  = e.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin( coneRadius_ ).first;
	iso2  = e.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin( coneRadius_ ).first;
	iso3  = e.isoDeposit(pat::PfGammaIso)->depositAndCountWithin( coneRadius_ ).first;
      }
      else {
	iso1   = e.dr03TkSumPt();//e.isoDeposit(tracker)->depositAndCountWithin(coneRadius_).first;//e.trackIso();
	iso2   = e.dr03EcalRecHitSumEt();
	iso3   = e.dr03HcalTowerSumEt();
      }
      
      iso4  = (iso1+iso2+iso3);
      //std::cout<<"ecalIso: "<<e.dr03EcalRecHitSumEt()<<std::endl;
      //std::cout<<"pt: "<<pt<<" TrackIso: "<<iso1<<", EcalIso: "<<iso2<<", HcalIso: "<<iso3<<std::endl;
 
     if(useRelIso_){
        //std::cout<<"singRelIso"<<std::endl;
	iso1 = iso1/pt;
	iso2 = iso2/pt;
	iso3 = iso3/pt;
	iso4 = iso4/pt;
      }
      //std::cout<<"pt: "<<pt<<" RelTrackIso: "<<iso1<<", RelEcalIso: "<<iso2<<", RelHcalIso: "<<iso3<<std::endl;
     
       if(e.isEB()){
	 cuts["ptx"]    = pt>ptX_;
	 cuts["etax"]   = fabs(eta)<etaX_;
	 cuts["sieie"]  = fabs(e.scSigmaIEtaIEta())<sieieBrlX_;
	 cuts["delphi"] = fabs(e.deltaPhiSuperClusterTrackAtVtx())<delphiBrlX_;
	 cuts["detain"] = fabs(e.deltaEtaSuperClusterTrackAtVtx())<detainBrlX_;
	 
 	 cuts["dcot"]   = fabs(e.convDcot())> dcotBrlX_;
         cuts["dist"]   = fabs(e.convDist())> distBrlX_;	
         cuts["hoe"]    = e.hadronicOverEm() <hoeBrlX_;
         cuts["nmshits"]= e.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits()<=nmisHitBrlX_;
	 
	 if(useCombIso_){
           cuts["iso1"]   = 1;
           cuts["iso2"]   = 1;
           cuts["iso3"]   = 1;
           cuts["iso4"]   = iso4 < iso4BrlX_;
         }
         else{
           cuts["iso1"]   = iso1 < iso1BrlX_;
           cuts["iso2"]   = iso2 < iso2BrlX_;
           cuts["iso3"]   = iso3 < iso3BrlX_;
           cuts["iso4"]   = 1;
           //std::cout<<"------------------------->"<<iso1<<"\t"<<iso2<<"\t"<<iso3<<std::endl;
           //std::cout<<"------------------------->"<<iso1BrlX_<<"\t"<<iso2BrlX_<<"\t"<<iso3BrlX_<<std::endl; 
           //std::cout<<"=========================>"<<cuts["iso1"]<<"\t"<<cuts["iso2"]<<"\t"<<cuts["iso3"]<<std::endl;

         }
	 
       }
       
       else if(e.isEE()){
	 cuts["ptx"]    = pt>ptX_;
	 cuts["etax"]   = fabs(e.eta())<etaX_;
	 cuts["sieie"]  = fabs(e.scSigmaIEtaIEta())<sieieEcpX_;
	 cuts["delphi"] = fabs(e.deltaPhiSuperClusterTrackAtVtx())<delphiEcpX_;
	 cuts["detain"] = fabs(e.deltaEtaSuperClusterTrackAtVtx())<detainEcpX_;
         cuts["dcot"]   = fabs(e.convDcot()) > dcotEcpX_;
         cuts["dist"]   = fabs(e.convDist()) > distEcpX_;
         cuts["hoe"]    = e.hadronicOverEm() < hoeEcpX_;
         cuts["nmshits"]= e.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits()<=nmisHitEcpX_;
	 if(useCombIso_){
	   cuts["iso1"]   = 1;
	   cuts["iso2"]   = 1;
	   cuts["iso3"]   = 1;
	   cuts["iso4"]   = iso4 < iso4EcpX_;
	 }
	 else{
	   cuts["iso1"]   = iso1 < iso1EcpX_;
	   cuts["iso2"]   = iso2 < iso2EcpX_;
	   cuts["iso3"]   = iso3 < iso3EcpX_;
	   cuts["iso4"]   = 1;
	 }
       }
      
    }
     
       else    {
 	 cuts["ptx"]    = 0;
	 cuts["etax"]   = 0;
	 cuts["sieie"]  = 0;
	 cuts["delphi"] = 0;
	 cuts["detain"] = 0;
  	 cuts["dcot"] = 0;
  	 cuts["dist"] = 0;	
         cuts["hoe"]    = 0;
         cuts["nmshits"] = 0;
	 cuts["iso1"]   = 0;
	 cuts["iso2"]   = 0;
	 cuts["iso3"]   = 0;
         cuts["iso4"]   = 0;
        
       }
   
       return cuts;
     }
