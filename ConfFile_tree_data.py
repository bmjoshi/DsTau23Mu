import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7')

process.load('RecoMET.METFilters.badGlobalMuonTaggersAOD_cff')
#switch on tagging mode:
process.badGlobalMuonTagger.taggingMode = cms.bool(True)
process.cloneGlobalMuonTagger.taggingMode = cms.bool(True)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/data/Run2016H/DoubleMuonLowMass/AOD/PromptReco-v2/000/281/613/00000/065CE7A9-E084-E611-B889-02163E012AAD.root',
    '/store/data/Run2016H/DoubleMuonLowMass/AOD/PromptReco-v2/000/281/613/00000/1AE06D85-EB84-E611-86A9-FA163E5F337F.root',
    #'/store/data/Run2016E/DoubleMuonLowMass/AOD/23Sep2016-v1/100000/F2C74B60-5494-E611-A55C-02163E0176BF.root'
    #'/store/data/Run2016B/DoubleMuonLowMass/AOD/23Sep2016-v3/00000/9856F893-AC9A-E611-88B4-FA163EECBDC5.root'
   ),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

)

from CondCore.DBCommon.CondDBSetup_cfi import *
process.l1Menu = cms.ESSource("PoolDBESSource",CondDBSetup,
                              connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                              toGet = cms.VPSet(cms.PSet(record = cms.string("L1TGlobalPrescalesVetosRcd"),
                                                         tag = cms.string("L1TGlobalPrescalesVetos_Stage2v0_hlt")),
                                                cms.PSet(record = cms.string("L1TUtmTriggerMenuRcd"),
                                                         tag = cms.string("L1TUtmTriggerMenu_Stage2v0_hlt"))
                                                )                              )
process.es_prefer_l1Menu = cms.ESPrefer("PoolDBESSource","l1Menu")


process.demo = cms.EDAnalyzer('dsTreeMaker',
    mid = cms.int32(15),
    MC = cms.bool(False),
    passhlt = cms.bool(True),
    wideSB = cms.bool(False),
    do2mu = cms.bool(True), # do 2 mu category or not
    muons = cms.InputTag("muons"),
    pvs = cms.InputTag("offlinePrimaryVertices"),
    svs = cms.InputTag("inclusiveSecondaryVertices"),
    trks = cms.InputTag("generalTracks"),
    btagsCvsB = cms.InputTag("pfCombinedCvsBJetTags"),
    btagsMVA = cms.InputTag("pfCombinedMVAV2BJetTags"),
    btagsCSV = cms.InputTag("pfCombinedSecondaryVertexV2BJetTags"),
    pfcands = cms.InputTag("particleFlow"),
    triggerBitsH = cms.InputTag("TriggerResults", "", "HLT"),
    triggerSummary = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    beamSpotHandle = cms.InputTag("offlineBeamSpot"),
    pileupSummary = cms.InputTag("addPileupInfo"),
    genParticles = cms.InputTag("genParticles"),
    AlgInputTag = cms.InputTag( "gtStage2Digis" ),
    BadGlbMuonFilter = cms.InputTag("badGlobalMuonTagger"),
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('tree_test.root'))

process.tagger = cms.Path(process.badGlobalMuonTagger)
process.p = cms.Path(process.demo)
process.schedule = cms.Schedule(process.tagger, process.p)

