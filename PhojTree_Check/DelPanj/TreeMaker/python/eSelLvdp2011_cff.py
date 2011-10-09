import FWCore.ParameterSet.Config as cms
from DelPanj.TreeMaker.eSelBase_cff import *

##The islation: non-combined, relative, detector based.
eIsoLvdpBrl = eIsoBase.clone()
eIsoLvdpBrl.iso1 = cms.double(0.09)##tkIso
eIsoLvdpBrl.iso2 = cms.double(0.07)##ecIso
eIsoLvdpBrl.iso3 = cms.double(0.10)##hcIso

eIsoLvdpEcp = eIsoLvdpBrl.clone()
eIsoLvdpEcp.iso1 = cms.double(0.04)##tkIso
eIsoLvdpEcp.iso2 = cms.double(0.05)##ecIso
eIsoLvdpEcp.iso3 = cms.double(0.025)##hcIso


##The are the gsf like id cuts, mva cut not included yet.
eIdLvdpBrl = eIdBase.clone()
eIdLvdpBrl.sieie = cms.double(0.01)
eIdLvdpBrl.delphi= cms.double(0.06)
eIdLvdpBrl.detain = cms.double(0.004)
eIdLvdpBrl.dist = cms.double(0.02)
eIdLvdpBrl.dcot = cms.double(0.02)
eIdLvdpBrl.hoe =  cms.double(0.04)
eIdLvdpBrl.nmisHit=cms.double(0.)

eIdLvdpEcp = eIdLvdpBrl.clone()
eIdLvdpEcp.sieie = cms.double(0.03)
eIdLvdpEcp.delphi= cms.double(0.03)
eIdLvdpEcp.detain = cms.double(0.007)
eIdLvdpEcp.hoe =  cms.double(0.15)


eSelLvdp2011 = eSelBase.clone()
eSelLvdp2011.coneRad= cms.double(0.3)
eSelLvdp2011.usePfIso= cms.bool(False)
eSelLvdp2011.useRelIso= cms.bool(True)
eSelLvdp2011.idBrl  = eIdLvdpBrl
eSelLvdp2011.idEcp  = eIdLvdpEcp
eSelLvdp2011.isoBrl = eIsoLvdpBrl
eSelLvdp2011.isoEcp = eIsoLvdpEcp
