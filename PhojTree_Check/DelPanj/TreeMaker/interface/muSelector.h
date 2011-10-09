#ifndef __MUSELECTOR_H__
#define __MUSELECTOR_H__

/*
MuonSelector
Optimal Usage: patMuons.

Lovedeep Kaur,
Panjab University.
*/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

class muSelector{
   
 public:
  muSelector(const edm::ParameterSet ps);
  std::map<std::string, bool> CutRecord(pat::Muon& mu, 
					reco::BeamSpot& beamSpot);
  
  ~muSelector(){}
  
  double ptX_;
  double etaX_;
  double dxyX_;
  double normalizedChi2X_;
  double trackerHitsX_;
  double pixelHitsX_;
  double muonHitsX_;
  double nMatchesX_;

  double iso1X_;
  double iso2X_;
  double iso3X_;
  double iso4X_;
  
  //some flags.
  bool usePfIso_;
  bool useRelIso_;
  bool useCombIso_;
  bool coneRadius_;
  bool applyGlbTrk_;    
   
  edm::ParameterSet idPar_;
  edm::ParameterSet isoPar_;
  


};

#endif
