#ifndef __ESELECTOR_H__
#define __ESELECTOR_H__

/*
ElectronSelector
Optimal Usage: pat Electrons.

Anil Singh,
Panjab University.
*/



#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/LorentzVector.h"

class eSelector{
   
 public:
  eSelector(const edm::ParameterSet ps);  
  std::map<std::string, bool> CutRecord(pat::Electron& e);
  
  
  
  ~eSelector(){}
  
  double ptX_;
  double etaX_;
  
  double sieieBrlX_;
  double delphiBrlX_;
  double detainBrlX_;
  double distBrlX_;
  double dcotBrlX_;
  double hoeBrlX_;
  double nmisHitBrlX_;
  double sieieEcpX_;
  double delphiEcpX_;
  double detainEcpX_;
  double distEcpX_;
  double dcotEcpX_;
  double hoeEcpX_;
  double nmisHitEcpX_;
  
  double iso1BrlX_;
  double iso2BrlX_;
  double iso3BrlX_;
  double iso4BrlX_;
  
  double iso1EcpX_;
  double iso2EcpX_;
  double iso3EcpX_;
  double iso4EcpX_;
  
  
   
  //some flags.
  bool usePfIso_;
  bool useRelIso_;
  bool useCombIso_;
  bool coneRadius_;
  
  
  //some pset declaration
  edm::ParameterSet idBrl_;
  edm::ParameterSet idEcp_;
  edm::ParameterSet isoBrl_;
  edm::ParameterSet isoEcp_; 
};

#endif
