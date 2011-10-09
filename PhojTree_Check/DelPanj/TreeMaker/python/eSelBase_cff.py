import FWCore.ParameterSet.Config as cms


###########################################################
#identification stuff                                     # 
###########################################################
eIdBase = cms.PSet(
    #Need to put in the variable values.
    sieie= cms.double(99.),
    delphi= cms.double(99.),
    detain=cms.double(99.),
    dist=cms.double(0.),
    dcot=cms.double(0.),
    hoe=cms.double(99.),
    nmisHit=cms.double(99.)
    )


##########################################################################
##usePfIso=True) iso1=chHadIso;iso2=neHadIso;iso3=gamIso;iso4=combIso    # 
##else iso1=tkIso;iso2=ecalIso;iso3=hcalIso;iso4=combIso                 #
##########################################################################
eIsoBase = cms.PSet(
    iso1 = cms.double(9999.),
    iso2 = cms.double(9999.),
    iso3 = cms.double(9999.),
    iso4 = cms.double(9999.)
    )


#######################################
#Base Electron Selector Configuration #
#######################################
eSelBase = cms.PSet(
    ptx= cms.double(10),
    etax= cms.double(2.3),
    ##Identification
    idBrl = eIdBase,
    idEcp = eIdBase,
    ##Isolation: Detector Based Rel Iso.
    coneRad= cms.double(0.3),
    usePfIso= cms.bool(False),
    useCombIso = cms.bool(False),
    useRelIso= cms.bool(True),
    isoBrl = eIsoBase,
    isoEcp = eIsoBase
    )
