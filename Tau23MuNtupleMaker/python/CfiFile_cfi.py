import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('Tau23MuNtupleMaker'
     ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
