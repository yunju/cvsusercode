import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:/data4/yunju/VBF2012/skim2l2q/SkimPAT_H2l2q_523_v3_l_GluGluToHToZZTo2L2Q_M-200h2l2qSkimData_22_1_2KD.root'
       #'file:/afs/cern.ch/work/y/yunju/public/SkimPAT_H2l2q_523_v3_l_GluGluToHToZZTo2L2Q_M-200h2l2qSkimData_22_1_2KD.root'
      'file:/afs/cern.ch/work/y/yunju/public/SkimPAT_H2l2q_533_v1_GluGluToHToZZTo2L2Q_M-200_h2l2qSkimData_2_1_5SE_2nd.root'
       

    )
)

process.demo = cms.EDAnalyzer('AnalyzerCheck'
)


process.p = cms.Path(process.demo)
