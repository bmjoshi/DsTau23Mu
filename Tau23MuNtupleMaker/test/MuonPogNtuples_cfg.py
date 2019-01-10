import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import subprocess
import sys

#from Filelist.filelist_2017_DoubleMuonLowMass_Run2017B_17Nov2017 import inputFiles 

options = VarParsing.VarParsing()

options.register('globalTag',
                 '94X_dataRun2_ReReco_EOY17_v6', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Global Tag")

options.register('nEvents',
                 -1, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Maximum number of processed events")

options.register('eosInputFolder',
                 'root://cmsdcadisk01.fnal.gov//dcache/uscmsdisk/store/relval/CMSSW_9_2_2/SingleMuon/RECO/91X_dataRun2_relval_v6_RelVal_sigMu2016E-v1/10000/806DF864-BF4D-E711-A60E-0025905B85E8.root',
                 #'/store/relval/CMSSW_9_2_0/RelValZMM_13/GEN-SIM-RECO/PU25ns_91X_mcRun2_asymptotic_v3-v1/10000', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "EOS folder with input files")

options.register('ntupleName',
                 'charmonium_test.root', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Folder and name ame for output ntuple")

options.register('runOnMC',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run on DATA or MC")

options.register('hltPathFilter',
                 "all", #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Filter on paths (now only accepts all or IsoMu20)")

options.register('minMuPt',
                 0.0, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Skim the ntuple saving only STA || TRK || GLB muons with pT > of this value")

options.register('minNMu',
                 2, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "number of TRK or GLB muons with pT > minMuPt to pass the skim")

options.parseArguments()

if options.hltPathFilter == "all" :
    pathCut   = "all"
    filterCut = "all"
elif options.hltPathFilter == "IsoMu20" :
    pathCut   = "HLT_IsoMu20_v"
    if options.runOnMC :
        filterCut = "hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09"
    else :
        filterCut = "hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09"
        
else :
    print "[" + sys.argv[0] + "]:", "hltPathFilter=", options.hltPathFilter, "is not a valid parameter!"
    sys.exit(100)

process = cms.Process("NTUPLES")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.nEvents))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = cms.string(options.globalTag)

process.source = cms.Source("PoolSource",                     
        fileNames = cms.untracked.vstring(),
        secondaryFileNames = cms.untracked.vstring()
)

process.output = cms.OutputModule("PoolOutputModule",
	     fileName = cms.untracked.string('tau23mu_2017b_crabTest.root'),
			)

#files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", options.eosInputFolder ])
#print inputFiles
#process.source.fileNames = [ options.eosInputFolder+"/"+f for f in files.split() ]  
#process.source.fileNames = inputFiles 
#process.source.fileNames = ['/store/data/Run2016D/DoubleMuonLowMass/AOD/07Aug17-v2/50000/8CFA50C6-259E-E711-AE20-00259075D72C.root','/store/data/Run2016D/DoubleMuonLowMass/AOD/07Aug17-v2/50000/A2C2EB1A-529F-E711-8FF7-0CC47AA53D68.root']
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

from MuonPOGtreeProducer.Tools.MuonPogNtuples_cff import appendMuonPogNtuple, customiseHlt, customiseMuonCuts
    
appendMuonPogNtuple(process,options.runOnMC,"HLT",options.ntupleName)
#customiseHlt(process,pathCut,filterCut)
#customiseMuonCuts(process,options.minMuPt,options.minNMu)
#process.prunedGenParticles.src = cms.InputTag("prunedGenParticles")
#process.goodOfflinePrimaryVertices.src = cms.InputTag("offlineSlimmedPrimaryVertices")
#process.MuonPogTree.GenTag = cms.untracked.InputTag("prunedGenParticles")
