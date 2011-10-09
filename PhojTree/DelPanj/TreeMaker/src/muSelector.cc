#include "DelPanj/TreeMaker/interface/muSelector.h"

muSelector::muSelector(const edm::ParameterSet ps):      
  ptX_  (ps.getParameter<double>("pt")),
  etaX_ (ps.getParameter<double>("eta")),

  idPar_(ps.getParameter<edm::ParameterSet> ("idPar")),
  isoPar_(ps.getParameter<edm::ParameterSet> ("isoPar")),

  //isolation parameters 
  usePfIso_(ps.getParameter<bool>("usePfIso")),  
  useRelIso_(ps.getParameter<bool>("useRelIso")),
  useCombIso_(ps.getParameter<bool>("useCombIso")),
  coneRadius_(ps.getParameter<double>("coneRad"))//,
{
  dxyX_ =(idPar_.getParameter<double>("dxy"));
  normalizedChi2X_=(idPar_.getParameter<double>("normalizedChi2"));
  trackerHitsX_=(idPar_.getParameter<double>("trackerHits"));
  pixelHitsX_=(idPar_.getParameter<double>("pixelHits"));
  muonHitsX_=(idPar_.getParameter<double>("muonHits"));
  nMatchesX_=(idPar_.getParameter<double>("nMatches"));
  applyGlbTrk_=(idPar_.getParameter<double>("applyGlbTrk"));//

  if(!useCombIso_){
    iso1X_ =(isoPar_.getParameter<double>("iso1"));
    iso2X_ =(isoPar_.getParameter<double>("iso2"));
    iso3X_ =(isoPar_.getParameter<double>("iso3"));
  }  
  else 
    iso4X_ =(isoPar_.getParameter<double>("iso4"));
}


std::map<std::string, bool>
muSelector::CutRecord(pat::Muon& mu, reco::BeamSpot& beamSpotHandle ){
  std::map<std::string, bool> cuts;
  double pt = mu.pt();
  double eta = mu.eta();
  double iso1 = 999.;
  double iso2 = 999.;
  double iso3 = 999.;
  double iso4 = 999.;
  
    if(usePfIso_){   
      iso1  = mu.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin( coneRadius_ ).first;
      iso2  = mu.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin( coneRadius_ ).first;
      iso3  = mu.isoDeposit(pat::PfGammaIso)->depositAndCountWithin( coneRadius_ ).first;
    }

    else {
      iso1   = mu.trackIso();
      iso2   = mu.ecalIso(); 
      iso3   = mu.hcalIso(); 
    }
   
    iso4  = (iso1+iso2+iso3);
    if(useRelIso_){
      iso1 = iso1/pt;
      iso2 = iso2/pt;
      iso3 = iso3/pt;
      iso4 = iso4/pt;
    }
    
    cuts["pt"]    = pt>ptX_;
    cuts["eta"]   = fabs(eta)<etaX_;
    if(applyGlbTrk_)
      cuts["isGlbNTrk"]   = mu.isGlobalMuon() && mu.isTrackerMuon(); 
    else
      cuts["isGlbNTrk"]   = 1;
    
    bool isTrackMuon =mu.isTrackerMuon();
    bool isGlobalMuon =mu.isGlobalMuon();
    bool isStandAloneMuon =mu.isStandAloneMuon();

    reco::TrackRef gm;
    if(isStandAloneMuon)
      gm = mu.standAloneMuon();
    if(isTrackMuon)
      gm = mu.track();
    if(isGlobalMuon)
      gm = mu.globalTrack();

    cuts["dxy"]           = fabs(gm->dxy(beamSpotHandle.position()))<dxyX_;
    cuts["normalizedChi2"]= gm->normalizedChi2()<normalizedChi2X_;
    cuts["trackerHits"]   = gm->hitPattern().numberOfValidTrackerHits()>trackerHitsX_; 
    cuts["pixelHits"]     = gm->hitPattern().numberOfValidPixelHits()>pixelHitsX_;
    cuts["muonHits"]      = gm->hitPattern().numberOfValidMuonHits()>muonHitsX_;	
    cuts["nMatches"]      = mu.numberOfMatches()>nMatchesX_;
    
    
    if(useCombIso_){
      cuts["iso1"]   = 1;
      cuts["iso2"]   = 1;
      cuts["iso3"]   = 1;
      cuts["iso4"]   = iso4 < iso4X_;
    }
    else{
      cuts["iso1"]   = iso1 < iso1X_;
      cuts["iso2"]   = iso2 < iso2X_;
      cuts["iso3"]   = iso3 < iso3X_;
      cuts["iso4"]   = 1;
    }
   
    return cuts;
}
