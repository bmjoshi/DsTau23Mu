// -*- C++ -*-
//
// Package:  DsTau3Mu/Tau3MuNtupleMaker/Tau23MuNtupleMaker.cc
// Class:  tau23mu
// 
// Original code: MuonPogNtuples.cc (Muon POG), dsTreeMaker.cc (Jian Wang)
// 

#include <iostream>
#include <algorithm>
#include <vector>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/src/classes.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "RecoMuon/MuonIdentification/interface/MuonCaloCompatibility.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/CTPPSDigi/interface/CTPPSPixelDigi.h"
#include "DataFormats/CTPPSDigi/interface/CTPPSPixelDigiCollection.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1TGlobalPrescalesVetosRcd.h"
#include "CondFormats/L1TObjects/interface/L1TGlobalPrescalesVetos.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetMatchInfo.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TTree.h"
#include "TLorentzVector.h"

#include "DsTau23Mu/Tau23MuNtupleMaker/src/Tau23MuTree.h"
#include "DsTau23Mu/Tau23MuNtupleMaker/src/Utils.h"
#include "DsTau23Mu/Tau23MuNtupleMaker/src/pdg_info.h"

//#include "RecoMuon/MuonIdentification/src/MuonKinkFinder.cc"
//#include "TrackingTools/TrackRefitter/interface/RefitDirection.h"
//#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
//#include "TrackingTools/TrackFitters/interface/TrajectoryFitterRecord.h"
//#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
//#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
//#include "scnew.h"

using namespace reco;
using namespace tau23mu;
using namespace l1t;
using namespace edm;
using namespace std;

#define dz_mumu_cut 0.5
#define dR_mumu_cut 0.8
#define K_mass_cut 0.03
#define triMuFitVtx_nC2_cut 5.0
#define diMuTrkFitVtx_nC2_cut 5.0
//#define min_fvnC_2mu1tk 10.0 

class Tau23MuNtupleMaker : public edm::EDAnalyzer 
{
    public:

      Tau23MuNtupleMaker(const edm::ParameterSet &);

      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void beginRun(const edm::Run&, const edm::EventSetup&);
      virtual void beginJob();
      virtual void endJob();

    private:
      void fillAnalysisTree(const edm::EventSetup&,
            const edm::Handle<edm::View<reco::Muon> > &,
            const edm::Handle<std::vector<reco::Track> >&,
            const edm::Handle<std::vector<reco::Vertex> > &,
            const edm::Handle<std::vector<reco::Vertex> > &,
            const edm::Handle<reco::BeamSpot> &,
            const edm::Handle<edm::TriggerResults> &, 
            const edm::Handle<trigger::TriggerEvent> &,
            const edm::TriggerNames &,
            bool ,bool);

      void fillGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > &,
            const  edm::Handle<GenEventInfoProduct> &);

      TransientVertex fitTriMuVertex(const edm::EventSetup&, TrackRef _trk1, TrackRef _trk2, TrackRef _trk3);

      void fillGenParticles(const edm::Handle<reco::GenParticleCollection> &);

      void fillHlt(const edm::Handle<edm::TriggerResults> &, 
            const edm::Handle<trigger::TriggerEvent> &,
            const edm::TriggerNames &);

      void fillPV(const edm::Handle<std::vector<reco::Vertex> > &);


      Int_t fillMuons(const edm::Handle<edm::View<reco::Muon> > &,
            const edm::Handle<std::vector<reco::Vertex> > &,
            const edm::Handle<reco::BeamSpot> &);

      void fillTracks(const edm::Handle<std::vector<reco::Track> >& ,
            const edm::Handle<std::vector<reco::Vertex> > &,
            const edm::Handle<reco::BeamSpot> &);

      void fillL1(const edm::Handle<l1t::MuonBxCollection> &);

      void printSummary();


      // returns false in case the match is for a RPC chamber
      bool getMuonChamberId(DetId & id, tau23mu::MuonDetType & det, Int_t & r, Int_t & phi, Int_t & eta) const ;

      edm::EDGetTokenT<edm::TriggerResults> trigResultsToken_;
      edm::EDGetTokenT<trigger::TriggerEvent> trigSummaryToken_;

      std::string trigFilterCut_;
      std::string trigPathCut_;

      edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
      edm::EDGetTokenT<std::vector<reco::Track> > trackToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > primaryVertexToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > secondaryVertexToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

      edm::EDGetTokenT<reco::PFMETCollection> pfMetToken_;
      edm::EDGetTokenT<reco::PFMETCollection> pfChMetToken_;
      edm::EDGetTokenT<reco::CaloMETCollection> caloMetToken_;

      edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileUpInfoToken_;
      edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;

      edm::EDGetTokenT<LumiScalersCollection> scalersToken_;

      edm::EDGetTokenT<l1t::MuonBxCollection> l1Token_;

      Float_t m_minMuPtCut;
      Int_t m_minNMuCut;

      tau23mu::Event event_;
      tau23mu::EventId eventId_;
      std::map<std::string,TTree*> tree_;

      Bool_t _runMC; 

      // branch flags
      Bool_t ev_doMuons;
      Bool_t ev_doTracks;
      Bool_t ev_doBS;
      Bool_t ev_doL1;
      Bool_t ev_doHLT;
      Bool_t ev_doVertexes;
      Bool_t ev_doMET;
      Bool_t ev_doTaus;
      Bool_t ev_doJets;
      Bool_t ev_doAnalysis;
      Bool_t ev_do2MuTrk;

};


Tau23MuNtupleMaker::Tau23MuNtupleMaker( const edm::ParameterSet & cfg )
{

    // Input collections
    edm::InputTag tag = cfg.getUntrackedParameter<edm::InputTag>("TrigResultsTag", edm::InputTag("TriggerResults::HLT"));
    if (tag.label() != "none") trigResultsToken_ = consumes<edm::TriggerResults>(tag);

    tag = cfg.getUntrackedParameter<edm::InputTag>("TrigSummaryTag", edm::InputTag("hltTriggerSummaryAOD::HLT")); 
    if (tag.label() != "none") trigSummaryToken_ =consumes<trigger::TriggerEvent>(tag);

    trigFilterCut_ = cfg.getUntrackedParameter<std::string>("TrigFilterCut", std::string("all"));
    trigPathCut_ = cfg.getUntrackedParameter<std::string>("TrigPathCut", std::string("all"));

    tag = cfg.getUntrackedParameter<edm::InputTag>("MuonTag", edm::InputTag("muons"));
    if (tag.label() != "none") muonToken_ = consumes<edm::View<reco::Muon> >(tag);

    tag = cfg.getUntrackedParameter<edm::InputTag>("PrimaryVertexTag", edm::InputTag("vertexes"));
    if (tag.label() != "none") primaryVertexToken_ = consumes<std::vector<reco::Vertex> >(tag);

    tag = cfg.getUntrackedParameter<edm::InputTag>("SecondaryVertexTag", edm::InputTag("secVertexes"));
    if (tag.label() != "none") secondaryVertexToken_ = consumes<std::vector<reco::Vertex> >(tag);

    tag = cfg.getUntrackedParameter<edm::InputTag>("TrackTag", edm::InputTag("tracks"));
    if (tag.label() != "none") trackToken_ = consumes<std::vector<reco::Track> >(tag);

    tag = cfg.getUntrackedParameter<edm::InputTag>("BeamSpotTag", edm::InputTag("offlineBeamSpot"));
    if (tag.label() != "none") beamSpotToken_ = consumes<reco::BeamSpot>(tag);

    tag = cfg.getUntrackedParameter<edm::InputTag>("PFMetTag", edm::InputTag("pfMet"));
    if (tag.label() != "none") pfMetToken_ = consumes<reco::PFMETCollection>(tag);

    tag = cfg.getUntrackedParameter<edm::InputTag>("PFChMetTag", edm::InputTag("pfChMet"));
    if (tag.label() != "none") pfChMetToken_ = consumes<reco::PFMETCollection>(tag);

    tag = cfg.getUntrackedParameter<edm::InputTag>("CaloMetTag", edm::InputTag("caloMet"));
    if (tag.label() != "none") caloMetToken_ = consumes<reco::CaloMETCollection>(tag); 

    tag = cfg.getUntrackedParameter<edm::InputTag>("GenTag", edm::InputTag("prunedGenParticles"));
    if (tag.label() != "none") genToken_ = consumes<reco::GenParticleCollection>(tag);

    tag = cfg.getUntrackedParameter<edm::InputTag>("PileUpInfoTag", edm::InputTag("pileupInfo"));
    if (tag.label() != "none") pileUpInfoToken_ = consumes<std::vector<PileupSummaryInfo> >(tag);

    tag = cfg.getUntrackedParameter<edm::InputTag>("GenInfoTag", edm::InputTag("generator"));
    if (tag.label() != "none") genInfoToken_ = consumes<GenEventInfoProduct>(tag);  

    tag = cfg.getUntrackedParameter<edm::InputTag>("ScalersTag", edm::InputTag("scalersRawToDigi"));
    if (tag.label() != "none") scalersToken_ = consumes<LumiScalersCollection>(tag);

    tag = cfg.getUntrackedParameter<edm::InputTag>("l1MuonsTag", edm::InputTag("gmtStage2Digis:Muon:"));
    if (tag.label() != "none") l1Token_ = consumes<l1t::MuonBxCollection>(tag);

    _runMC = cfg.getUntrackedParameter<bool>("mc_",false);

    m_minMuPtCut = cfg.getUntrackedParameter<double>("MinMuPtCut", 0.);
    m_minNMuCut  = cfg.getUntrackedParameter<int>("MinNMuCut",  0.);

    ev_doMuons = cfg.getUntrackedParameter<bool>("_doMuons", false);
    ev_doTracks = cfg.getUntrackedParameter<bool>("_doTracks", false);
    ev_doBS = cfg.getUntrackedParameter<bool>("_doBS", false);
    ev_doL1 = cfg.getUntrackedParameter<bool>("_doL1", false);
    ev_doHLT = cfg.getUntrackedParameter<bool>("_doHLT", false);
    ev_doVertexes = cfg.getUntrackedParameter<bool>("_doVertexes", false);
    ev_doMET = cfg.getUntrackedParameter<bool>("_doMET", false);
    ev_doJets = cfg.getUntrackedParameter<bool>("_doJets", false);
    ev_doAnalysis = cfg.getUntrackedParameter<bool>("_doAnalysis",true);
    ev_doTaus = cfg.getUntrackedParameter<bool>("_doTaus", false);
    ev_do2MuTrk = cfg.getUntrackedParameter<bool>("_do2MuTrk",true);

}


void Tau23MuNtupleMaker::beginJob() 
{

    edm::Service<TFileService> fs;
    tree_["tau23muTree"] = fs->make<TTree>("TAU23MUTREE","Tau to 3 Mu tree");

    int splitBranches = 2;
    tree_["tau23muTree"]->Branch("event",&event_,64000,splitBranches);
    tree_["tau23muTree"]->Branch("eventId",&eventId_,64000,splitBranches);

}


void Tau23MuNtupleMaker::beginRun(const edm::Run & run, const edm::EventSetup & config )
{

}


void Tau23MuNtupleMaker::endJob() 
{

}


void Tau23MuNtupleMaker::analyze (const edm::Event & ev, const edm::EventSetup & iSetup)
{

    // Clearing branch variables
    // and setting default values
    event_.hlt.triggers.clear();
    event_.hlt.objects.clear();
    event_.l1muons.clear();

    event_.genParticles.clear();
    event_.genInfos.clear();
    event_.muons.clear();
    event_.tracks.clear();
    event_.triMuon_cand_coll.clear();
    event_.diMuonTrk_cand_coll.clear();

    event_.mets.pfMet = -999; 
    event_.mets.pfChMet = -999; 
    event_.mets.caloMet = -999; 

    for (unsigned int ix=0; ix<3; ++ix) {
      event_.primaryVertex[ix] = 0.;
      for (unsigned int iy=0; iy<3; ++iy) {
        event_.cov_primaryVertex[ix][iy] = 0.;
      }
    }

    event_.nVtx = -1;
    event_.nTrk = -1;
    event_.n3mu = -1;
    event_.n2muTrk = -1;

    // Fill general information
    // run, luminosity block, event
    event_.runNumber = ev.id().run();
    event_.luminosityBlockNumber = ev.id().luminosityBlock();
    event_.eventNumber = ev.id().event();

    if (_runMC) edm::LogInfo("") <<"[Tau23MuNtupleMaker]: Running on CMS dataset...";
    else  edm::LogInfo("") <<"[Tau23MuNtupleMaker]: Running on MC..."; 

    // Fill GEN pile up information
    edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
    edm::Handle<GenEventInfoProduct> genInfo;
    if (!ev.isRealData()) 
    {
      if (!pileUpInfoToken_.isUninitialized() &&
            !genInfoToken_.isUninitialized()) 
      {
        if (ev.getByToken(pileUpInfoToken_, puInfo) &&
              ev.getByToken(genInfoToken_, genInfo) ) 
            fillGenInfo(puInfo,genInfo);
        else 
            edm::LogError("") << "[Tau23MuNtupleMaker]: Pile-Up Info collection does not exist !!!";
      }  
    }


    // Fill GEN particles information
    edm::Handle<reco::GenParticleCollection> genParticles;
    if (!ev.isRealData()) 
    {
      if (!genToken_.isUninitialized() ) 
      { 
        if (ev.getByToken(genToken_, genParticles)) 
            fillGenParticles(genParticles);
        else 
            edm::LogError("") << "[Tau23MuNtupleMaker]: GEN collection does not exist !!!";
      }
    }

    edm::Handle<LumiScalersCollection> lumiScalers;
    if (ev.isRealData()) 
    {
      event_.bxId  = ev.bunchCrossing();
      event_.orbit = ev.orbitNumber();
      if (!scalersToken_.isUninitialized()) 
      { 
        if (ev.getByToken(scalersToken_, lumiScalers) && 
              lumiScalers->size() > 0 ) 
            event_.instLumi  = lumiScalers->begin()->instantLumi();
        else 
            edm::LogError("") << "[Tau23MuNtupleMaker]: Scaler collection does not exist !!!";
      }
    }
    edm::Handle<edm::TriggerResults> triggerResults;
    edm::Handle<trigger::TriggerEvent> triggerEvent;
    // Fill trigger information
    if (!trigResultsToken_.isUninitialized() &&
        !trigSummaryToken_.isUninitialized()) 
    {
      if (ev.getByToken(trigResultsToken_, triggerResults) &&
            ev.getByToken(trigSummaryToken_, triggerEvent)) 
        fillHlt(triggerResults, triggerEvent,ev.triggerNames(*triggerResults));
      else 
        edm::LogError("") << "[Tau23MuNtupleMaker]: Trigger collections do not exist !!!";
    }


    // Fill vertex information
    edm::Handle<std::vector<reco::Vertex> > vertexes;
    if(!primaryVertexToken_.isUninitialized()) 
    {
      if (ev.getByToken(primaryVertexToken_, vertexes))
        fillPV(vertexes);
      else 
        edm::LogError("") << "[Tau23MuNtupleMaker]: Vertex collection does not exist !!!";
    }

    // Fill secondary vertex information
    edm::Handle<std::vector<reco::Vertex> > secVertexes;
    if(!primaryVertexToken_.isUninitialized()) 
    {
      if (ev.getByToken(secondaryVertexToken_, secVertexes))
        fillPV(secVertexes);
      else 
        edm::LogError("") << "[Tau23MuNtupleMaker]: Secondary vertex collection does not exist !!!";
    }


    // Get beam spot for muons
    edm::Handle<reco::BeamSpot> beamSpot;
    if (!beamSpotToken_.isUninitialized() && ev_doMET) 
    { 
      if (!ev.getByToken(beamSpotToken_, beamSpot)) 
        edm::LogError("") << "[Tau23MuNtupleMaker]: Beam spot collection not found !!!";
    }

    // Fill (raw) MET information: PF, PF charged, Calo  
    edm::Handle<reco::PFMETCollection> pfMet; 
    if(!pfMetToken_.isUninitialized() && ev_doMET) 
    { 
      if (!ev.getByToken(pfMetToken_, pfMet)) 
        edm::LogError("") << "[Tau23MuNtupleMaker]: PFMet collection does not exist !!!"; 
      else { 
        const reco::PFMET &iPfMet = (*pfMet)[0]; 
        event_.mets.pfMet = iPfMet.et(); 
      } 
    } 

    edm::Handle<reco::PFMETCollection> pfChMet; 
    if(!pfChMetToken_.isUninitialized() && ev_doMET) 
    { 
      if (!ev.getByToken(pfChMetToken_, pfChMet)) 
        edm::LogError("") << "[Tau23MuNtupleMaker]: PFChMet collection does not exist !!!"; 
      else { 
        const reco::PFMET &iPfChMet = (*pfChMet)[0]; 
        event_.mets.pfChMet = iPfChMet.et(); 
      } 
    } 

    edm::Handle<reco::CaloMETCollection> caloMet; 
    if(!caloMetToken_.isUninitialized() && ev_doMET) 
    { 
      if (!ev.getByToken(caloMetToken_, caloMet)) 
        edm::LogError("") << "[Tau23MuNtupleMaker]: CaloMet collection does not exist !!!"; 
      else { 
        const reco::CaloMET &iCaloMet = (*caloMet)[0]; 
        event_.mets.caloMet = iCaloMet.et(); 
      } 
    } 

    // Get muons  
    edm::Handle<edm::View<reco::Muon> > muons;
    if (!muonToken_.isUninitialized() && ev_doMuons) 
    { 
      if (!ev.getByToken(muonToken_, muons)) 
        edm::LogError("") << "[Tau23MuNtupleMaker]: Muon collection does not exist !!!";
    }

    //Get tracks
    edm::Handle<std::vector<reco::Track> > tracks;
    if (!trackToken_.isUninitialized() && ev_doTracks)
    {
      if (!ev.getByToken(trackToken_, tracks))
        edm::LogError("") << "[Tau23MuNtupleMaker]: Track collection does not exist !!!";
      else fillTracks(tracks,vertexes,beamSpot);

    }
    //Get Jet ??
    //Get secondary vertices ??

    Int_t nGoodMuons = 0;
    eventId_.maxPTs.clear();
    // Fill muon information
    if (muons.isValid() && vertexes.isValid() && beamSpot.isValid() && ev_doMuons) 
    {
      nGoodMuons = fillMuons(muons,vertexes,beamSpot);
    }
    eventId_.nMuons = nGoodMuons;


    //Fill L1 informations
    edm::Handle<l1t::MuonBxCollection> l1s;
    if (!l1Token_.isUninitialized() && ev_doL1)
    {
      if (!ev.getByToken(l1Token_, l1s))
        edm::LogError("") << "[Tau23MuNtupleMaker] L1 muon bx collection does not exist !!!";
      else {
        fillL1(l1s);
      }
    }
    // Fill 3mu candidates and 2mu + trk candidates
    if (muons.isValid() && vertexes.isValid() && beamSpot.isValid() && secVertexes.isValid() && tracks.isValid() && ev_doAnalysis)
    {
      if (ev.getByToken(l1Token_, l1s) && ev.getByToken(trackToken_, tracks) 
            && ev.getByToken(muonToken_, muons) && ev.getByToken(beamSpotToken_, beamSpot) &&
            ev.getByToken(primaryVertexToken_, vertexes) && ev.getByToken(secondaryVertexToken_, secVertexes) &&
            ev.getByToken(trigResultsToken_, triggerResults) && ev.getByToken(trigSummaryToken_, triggerEvent)){

        fillAnalysisTree(iSetup, muons,tracks,vertexes,secVertexes,beamSpot,triggerResults,triggerEvent,ev.triggerNames(*triggerResults),_runMC,ev_do2MuTrk);

      }

      else edm::LogError("") <<"[Tau23MuNtupleMaker]: One of the collections needed for analysis does not exist !!!";
    }



    //if (nGoodMuons >= m_minNMuCut)
    tree_["tau23muTree"]->Fill();

}

void Tau23MuNtupleMaker::fillGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > & puInfo,
      const edm::Handle<GenEventInfoProduct> & gen)
{

    tau23mu::GenInfo genInfo;

    genInfo.trueNumberOfInteractions = -1.;
    genInfo.actualNumberOfInteractions = -1.;
    genInfo.genWeight = gen->weight() ;

    std::vector<PileupSummaryInfo>::const_iterator puInfoIt  = puInfo->begin();
    std::vector<PileupSummaryInfo>::const_iterator puInfoEnd = puInfo->end();

    for(; puInfoIt != puInfoEnd; ++puInfoIt) 
    {
      int bx = puInfoIt->getBunchCrossing();

      if(bx == 0) 
      { 
        genInfo.trueNumberOfInteractions = puInfoIt->getTrueNumInteractions();
        genInfo.actualNumberOfInteractions = puInfoIt->getPU_NumInteractions();
        continue;
      }
    }

    event_.genInfos.push_back(genInfo);

}


void Tau23MuNtupleMaker::fillGenParticles(const edm::Handle<reco::GenParticleCollection> & genParticles)
{

    unsigned int gensize = genParticles->size();

    // Do not record the initial protons
    for (unsigned int i=0; i<gensize; ++i) 
    {

      const reco::GenParticle& part = genParticles->at(i);

      tau23mu::GenParticle gensel;
      gensel.pdgId = part.pdgId();
      gensel.status = part.status();
      gensel.energy = part.energy();
      gensel.pt = part.pt();
      gensel.eta = part.eta();
      gensel.phi = part.phi();
      gensel.vx = part.vx();
      gensel.vy = part.vy();
      gensel.vz = part.vz();

      // Full set of GenFlags
      gensel.flags.clear();
      reco::GenStatusFlags statusflags = part.statusFlags();
      if (statusflags.flags_.size() == 15)
        for (unsigned int flag = 0; flag < statusflags.flags_.size(); ++flag)
            gensel.flags.push_back(statusflags.flags_[flag]);  

      gensel.mothers.clear();
      unsigned int nMothers = part.numberOfMothers();

      for (unsigned int iMother=0; iMother<nMothers; ++iMother) 
      {
        gensel.mothers.push_back(part.motherRef(iMother)->pdgId());
      }

      gensel.daughters.clear();
      unsigned int nDaughters = part.numberOfDaughters();

      for (unsigned int iDaughter=0; iDaughter<nDaughters; ++iDaughter)
      {
        gensel.daughters.push_back(part.daughterRef(iDaughter)->pdgId());
      }

      // Protect agains bug in genParticles (missing mother => first proton)
      if (i>=2 && nMothers==0) gensel.mothers.push_back(0);

      //Protects against childless genParticles
      if (nDaughters==0) gensel.daughters.push_back(0);

      event_.genParticles.push_back(gensel);
    }

}


void Tau23MuNtupleMaker::fillHlt(const edm::Handle<edm::TriggerResults> & triggerResults, 
      const edm::Handle<trigger::TriggerEvent> & triggerEvent,
      const edm::TriggerNames & triggerNames)
{  

    for (unsigned int iTrig=0; iTrig<triggerNames.size(); ++iTrig) 
    {

      if (triggerResults->accept(iTrig)) 
      {
        std::string pathName = triggerNames.triggerName(iTrig);
        if (trigPathCut_ == "all" || pathName.find(trigPathCut_) != std::string::npos)
            event_.hlt.triggers.push_back(pathName);
      }
    }

    const trigger::size_type nFilters(triggerEvent->sizeFilters());

    for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter) 
    {

      std::string filterTag = triggerEvent->filterTag(iFilter).encode();

      if (trigFilterCut_ == "all" || filterTag.find(trigFilterCut_) != std::string::npos)
      {

        trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
        const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());

        for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey) 
        {  
            trigger::size_type objKey = objectKeys.at(iKey);
            const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);

            tau23mu::HLTObject hltObj;

            float trigObjPt = triggerObj.pt();
            float trigObjEta = triggerObj.eta();
            float trigObjPhi = triggerObj.phi();

            hltObj.filterTag = filterTag;

            hltObj.pt  = trigObjPt;
            hltObj.eta = trigObjEta;
            hltObj.phi = trigObjPhi;

            event_.hlt.objects.push_back(hltObj);

        }
      }
    }

}

void Tau23MuNtupleMaker::fillL1(const edm::Handle<l1t::MuonBxCollection> & l1MuonBxColl)
{

    for (int ibx = l1MuonBxColl->getFirstBX(); ibx <= l1MuonBxColl->getLastBX(); ++ibx) 
    {
      for (auto l1MuIt = l1MuonBxColl->begin(ibx); l1MuIt != l1MuonBxColl->end(ibx); ++l1MuIt)
      {

        tau23mu::L1Muon l1part;
        l1part.pt = l1MuIt->pt();
        l1part.eta = l1MuIt->eta();
        l1part.phi = l1MuIt->phi();
        l1part.charge = l1MuIt->hwChargeValid() ? l1MuIt->charge() : 0;

        l1part.quality = l1MuIt->hwQual();
        l1part.bx = ibx;

        l1part.tfIndex = l1MuIt->tfMuonIndex();

        event_.l1muons.push_back(l1part);

      }
    }
}



void Tau23MuNtupleMaker::fillPV(const edm::Handle<std::vector<reco::Vertex> > & vertexes)
{

    int nVtx = 0;

    std::vector<reco::Vertex>::const_iterator vertexIt  = vertexes->begin();
    std::vector<reco::Vertex>::const_iterator vertexEnd = vertexes->end();

    for (; vertexIt != vertexEnd; ++vertexIt) 
    {

      const reco::Vertex& vertex = *vertexIt;

      if (!vertex.isValid()) continue;
      ++nVtx;

      if (vertexIt == vertexes->begin()) 
      {
        event_.primaryVertex[0] = vertex.x();
        event_.primaryVertex[1] = vertex.y();
        event_.primaryVertex[2] = vertex.z();

        for (unsigned int ix=0; ix<3; ++ix) 
        {
            for (unsigned int iy=0; iy<3; ++iy) 
            {
              event_.cov_primaryVertex[ix][iy] = vertex.covariance(ix,iy);
            }
        }
      }
    }

    event_.nVtx = nVtx;
}

Int_t Tau23MuNtupleMaker::fillMuons(const edm::Handle<edm::View<reco::Muon> > & muons,
      const edm::Handle<std::vector<reco::Vertex> > & vertexes,
      const edm::Handle<reco::BeamSpot> & beamSpot){

    edm::View<reco::Muon>::const_iterator muonIt  = muons->begin();
    edm::View<reco::Muon>::const_iterator muonEnd = muons->end();

    for (; muonIt != muonEnd; ++muonIt) {

      const reco::Muon& mu = (*muonIt);

      bool isGlobal  = mu.isGlobalMuon();
      bool isTracker = mu.isTrackerMuon();
      bool isTrackerArb  = muon::isGoodMuon(mu, muon::TrackerMuonArbitrated); 
      bool isRPC = mu.isRPCMuon();
      bool isStandAlone  = mu.isStandAloneMuon();
      bool isPF = mu.isPFMuon();

      bool hasInnerTrack = !mu.innerTrack().isNull();
      bool hasTunePTrack = !mu.tunePMuonBestTrack().isNull();
      bool hasPickyTrack = !mu.pickyTrack().isNull();
      bool hasDytTrack = !mu.dytTrack().isNull();
      bool hasTpfmsTrack = !mu.tpfmsTrack().isNull();

      tau23mu::Muon ntupleMu;

      ntupleMu.pt = mu.pt();
      ntupleMu.eta  = mu.eta();
      ntupleMu.phi  = mu.phi();
      ntupleMu.charge = mu.charge();
      ntupleMu.vx = mu.vx();
      ntupleMu.vy = mu.vy();
      ntupleMu.vz = mu.vz();

      ntupleMu.fits.push_back(tau23mu::MuonFit(mu.pt(),mu.eta(),mu.phi(),
              mu.charge(),mu.muonBestTrack()->ptError()));

      ntupleMu.fits.push_back(tau23mu::MuonFit(hasInnerTrack ? mu.innerTrack()->pt()  : -1000.,
              hasInnerTrack ? mu.innerTrack()->eta() : -1000.,
              hasInnerTrack ? mu.innerTrack()->phi() : -1000.,
              hasInnerTrack ? mu.innerTrack()->charge()  : -1000.,
              hasInnerTrack ? mu.innerTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(tau23mu::MuonFit(isStandAlone ? mu.outerTrack()->pt()  : -1000.,
              isStandAlone ? mu.outerTrack()->eta() : -1000.,
              isStandAlone ? mu.outerTrack()->phi() : -1000.,
              isStandAlone ? mu.outerTrack()->charge()  : -1000.,
              isStandAlone ? mu.outerTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(tau23mu::MuonFit(isGlobal ? mu.globalTrack()->pt()  : -1000.,
              isGlobal ? mu.globalTrack()->eta() : -1000.,
              isGlobal ? mu.globalTrack()->phi() : -1000.,
              isGlobal ? mu.globalTrack()->charge()  : -1000.,
              isGlobal ? mu.globalTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(tau23mu::MuonFit(hasTunePTrack ? mu.tunePMuonBestTrack()->pt()  : -1000.,
              hasTunePTrack ? mu.tunePMuonBestTrack()->eta() : -1000.,
              hasTunePTrack ? mu.tunePMuonBestTrack()->phi() : -1000.,
              hasTunePTrack ? mu.tunePMuonBestTrack()->charge()  : -1000.,
              hasTunePTrack ? mu.tunePMuonBestTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(tau23mu::MuonFit(hasPickyTrack ? mu.pickyTrack()->pt()  : -1000.,
              hasPickyTrack ? mu.pickyTrack()->eta() : -1000.,
              hasPickyTrack ? mu.pickyTrack()->phi() : -1000.,
              hasPickyTrack ? mu.pickyTrack()->charge()  : -1000.,
              hasPickyTrack ? mu.pickyTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(tau23mu::MuonFit(hasDytTrack ? mu.dytTrack()->pt()  : -1000.,
              hasDytTrack ? mu.dytTrack()->eta() : -1000.,
              hasDytTrack ? mu.dytTrack()->phi() : -1000.,
              hasDytTrack ? mu.dytTrack()->charge()  : -1000.,
              hasDytTrack ? mu.dytTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(tau23mu::MuonFit(hasTpfmsTrack ? mu.tpfmsTrack()->pt()  : -1000.,
              hasTpfmsTrack ? mu.tpfmsTrack()->eta() : -1000.,
              hasTpfmsTrack ? mu.tpfmsTrack()->phi() : -1000.,
              hasTpfmsTrack ? mu.tpfmsTrack()->charge()  : -1000.,
              hasTpfmsTrack ? mu.tpfmsTrack()->ptError() : -1000.));

      // Detector Based Isolation
      reco::MuonIsolation detIso03 = mu.isolationR03();

      ntupleMu.trackerIso = detIso03.sumPt;
      ntupleMu.EMCalIso = detIso03.emEt;
      ntupleMu.HCalIso  = detIso03.hadEt;

      // PF Isolation
      reco::MuonPFIsolation pfIso04 = mu.pfIsolationR04();
      reco::MuonPFIsolation pfIso03 = mu.pfIsolationR03();

      ntupleMu.chargedHadronIso = pfIso04.sumChargedHadronPt;
      ntupleMu.chargedHadronIsoPU = pfIso04.sumPUPt; 
      ntupleMu.neutralHadronIso = pfIso04.sumNeutralHadronEt;
      ntupleMu.photonIso = pfIso04.sumPhotonEt;

      ntupleMu.isGlobal = isGlobal ? 1 : 0; 
      ntupleMu.isTracker  = isTracker ? 1 : 0; 
      ntupleMu.isTrackerArb = isTrackerArb ? 1 : 0; 
      ntupleMu.isRPC = isRPC ? 1 : 0;
      ntupleMu.isStandAlone = isStandAlone ? 1 : 0;
      ntupleMu.isPF = isPF ? 1 : 0;

      ntupleMu.nHitsGlobal = isGlobal ? mu.globalTrack()->numberOfValidHits() : -999; 
      ntupleMu.nHitsTracker  = isTracker  ? mu.innerTrack()->numberOfValidHits()  : -999; 
      ntupleMu.nHitsStandAlone = isStandAlone ? mu.outerTrack()->numberOfValidHits()  : -999;

      ntupleMu.glbNormChi2 = isGlobal  ? mu.globalTrack()->normalizedChi2() : -999; 
      ntupleMu.trkNormChi2 = hasInnerTrack ? mu.innerTrack()->normalizedChi2()  : -999; 
      ntupleMu.trkMuonMatchedStations = isTracker ? mu.numberOfMatchedStations() : -999; 
      ntupleMu.glbMuonValidHits = isGlobal  ? mu.globalTrack()->hitPattern().numberOfValidMuonHits() : -999; 
      ntupleMu.trkPixelValidHits = hasInnerTrack ? mu.innerTrack()->hitPattern().numberOfValidPixelHits() : -999; 
      ntupleMu.trkPixelLayersWithMeas = hasInnerTrack ? mu.innerTrack()->hitPattern().pixelLayersWithMeasurement() : -999; 
      ntupleMu.trkTrackerLayersWithMeas = hasInnerTrack ? mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() : -999; 

      ntupleMu.bestMuPtErr = mu.muonBestTrack()->ptError(); 

      ntupleMu.trkValidHitFrac = hasInnerTrack  ? mu.innerTrack()->validFraction() : -999; 
      ntupleMu.trkStaChi2  = isGlobal ? mu.combinedQuality().chi2LocalPosition : -999; 
      ntupleMu.trkKink = isGlobal ? mu.combinedQuality().trkKink  : -999; 
      ntupleMu.muSegmComp  = (isGlobal || isTracker) ? muon::segmentCompatibility(mu) : -999; 

      ntupleMu.isTrkMuOST  = muon::isGoodMuon(mu, muon::TMOneStationTight) ? 1 : 0; 
      ntupleMu.isTrkHP = hasInnerTrack && mu.innerTrack()->quality(reco::TrackBase::highPurity) ? 1 : 0; 

      if ( mu.isMatchesValid() && ntupleMu.isTrackerArb )
      {
        for ( reco::MuonChamberMatch match : mu.matches() )
        {
            tau23mu::ChambMatch ntupleMatch;

            if ( getMuonChamberId(match.id,
                    ntupleMatch.type,ntupleMatch.r,
                    ntupleMatch.phi,ntupleMatch.eta)
             )
            {

              ntupleMatch.x = mu.trackX(match.station(),match.detector());
              ntupleMatch.y = mu.trackY(match.station(),match.detector());
              ntupleMatch.dXdZ = mu.trackDxDz(match.station(),match.detector());
              ntupleMatch.dYdZ = mu.trackDyDz(match.station(),match.detector());

              ntupleMatch.errxTk = mu.trackXErr(match.station(),match.detector());
              ntupleMatch.erryTk = mu.trackYErr(match.station(),match.detector());

              ntupleMatch.errDxDzTk = mu.trackDxDzErr(match.station(),match.detector());
              ntupleMatch.errDyDzTk = mu.trackDyDzErr(match.station(),match.detector());

              ntupleMatch.dx = mu.dX(match.station(),match.detector());
              ntupleMatch.dy = mu.dY(match.station(),match.detector());
              ntupleMatch.dDxDz = mu.dDxDz(match.station(),match.detector());
              ntupleMatch.dDyDz = mu.dDxDz(match.station(),match.detector());

              ntupleMatch.errxSeg = mu.segmentXErr(match.station(),match.detector());
              ntupleMatch.errySeg = mu.segmentYErr(match.station(),match.detector());
              ntupleMatch.errDxDzSeg = mu.segmentDxDzErr(match.station(),match.detector());
              ntupleMatch.errDyDzSeg = mu.segmentDyDzErr(match.station(),match.detector());

              ntupleMu.matches.push_back(ntupleMatch);
            }
        }
      }

      ntupleMu.dxyBest  = -999; 
      ntupleMu.dzBest = -999; 
      ntupleMu.dxyInner = -999; 
      ntupleMu.dzInner  = -999; 

      ntupleMu.isoPflow04 = (pfIso04.sumChargedHadronPt+ 
            std::max(0.,pfIso04.sumPhotonEt+pfIso04.sumNeutralHadronEt - 0.5*pfIso04.sumPUPt)) / mu.pt();

      ntupleMu.isoPflow03 = (pfIso03.sumChargedHadronPt+ 
            std::max(0.,pfIso03.sumPhotonEt+pfIso03.sumNeutralHadronEt - 0.5*pfIso03.sumPUPt)) / mu.pt();

      double dxybs = isGlobal ? mu.globalTrack()->dxy(beamSpot->position()) :
        hasInnerTrack ? mu.innerTrack()->dxy(beamSpot->position()) : -1000;
      double dzbs  = isGlobal ? mu.globalTrack()->dz(beamSpot->position()) :
        hasInnerTrack ? mu.innerTrack()->dz(beamSpot->position()) : -1000;

      double dxy = -1000.;
      double dz  = -1000.;

      ntupleMu.isSoft  = 0; 
      ntupleMu.isTight = 0; 
      ntupleMu.isHighPt  = 0;
      ntupleMu.isLoose = muon::isLooseMuon(mu)  ? 1 : 0; 
      ntupleMu.isMedium  = muon::isMediumMuon(mu) ? 1 : 0; 

      if (vertexes->size() > 0)
      {
        const reco::Vertex & vertex = vertexes->at(0);

        dxy = isGlobal ? mu.globalTrack()->dxy(vertex.position()) :
            hasInnerTrack ? mu.innerTrack()->dxy(vertex.position()) : -1000;
        dz = isGlobal ? mu.globalTrack()->dz(vertex.position()) :
            hasInnerTrack ? mu.innerTrack()->dz(vertex.position()) : -1000;

        ntupleMu.dxyBest  = mu.muonBestTrack()->dxy(vertex.position()); 
        ntupleMu.dzBest = mu.muonBestTrack()->dz(vertex.position()); 
        if(hasInnerTrack) { 
            ntupleMu.dxyInner = mu.innerTrack()->dxy(vertex.position()); 
            ntupleMu.dzInner  = mu.innerTrack()->dz(vertex.position()); 
        } 

        ntupleMu.isSoft  = muon::isSoftMuon(mu,vertex) ? 1 : 0; 
        ntupleMu.isTight = muon::isTightMuon(mu,vertex)  ? 1 : 0; 
        ntupleMu.isHighPt  = muon::isHighPtMuon(mu,vertex) ? 1 : 0;

      }

      ntupleMu.dxy  = dxy;
      ntupleMu.dz = dz;
      ntupleMu.edxy = isGlobal ? mu.globalTrack()->dxyError() : hasInnerTrack ? mu.innerTrack()->dxyError() : -1000;
      ntupleMu.edz  = isGlobal ? mu.globalTrack()->dzError()  : hasInnerTrack ? mu.innerTrack()->dzError() : -1000;

      ntupleMu.dxybs  = dxybs;
      ntupleMu.dzbs = dzbs;

      if(mu.isTimeValid()) { 
        ntupleMu.muonTimeDof = mu.time().nDof; 
        ntupleMu.muonTime  = mu.time().timeAtIpInOut; 
        ntupleMu.muonTimeErr = mu.time().timeAtIpInOutErr; 
      } 
      else { 
        ntupleMu.muonTimeDof = -999; 
        ntupleMu.muonTime  = -999; 
        ntupleMu.muonTimeErr = -999; 
      } 

      if(mu.rpcTime().nDof > 0) { 
        ntupleMu.muonRpcTimeDof = mu.rpcTime().nDof; 
        ntupleMu.muonRpcTime  = mu.rpcTime().timeAtIpInOut; 
        ntupleMu.muonRpcTimeErr = mu.rpcTime().timeAtIpInOutErr; 
      } 
      else { 
        ntupleMu.muonRpcTimeDof = -999; 
        ntupleMu.muonRpcTime  = -999; 
        ntupleMu.muonRpcTimeErr = -999; 
      } 

      // asking for a TRK or GLB muon with minimal pT cut
      // ignoring STA muons in this logic
      if ( m_minMuPtCut < 0 ||
            (
           (isTracker || isGlobal || isStandAlone) &&
           (ntupleMu.fitPt(tau23mu::MuonFitType::DEFAULT) > m_minMuPtCut ||
            ntupleMu.fitPt(tau23mu::MuonFitType::GLB) > m_minMuPtCut ||
            ntupleMu.fitPt(tau23mu::MuonFitType::TUNEP) > m_minMuPtCut ||
            ntupleMu.fitPt(tau23mu::MuonFitType::INNER) > m_minMuPtCut ||
            ntupleMu.fitPt(tau23mu::MuonFitType::STA) > m_minMuPtCut ||
            ntupleMu.fitPt(tau23mu::MuonFitType::PICKY) > m_minMuPtCut ||
            ntupleMu.fitPt(tau23mu::MuonFitType::DYT) > m_minMuPtCut ||
            ntupleMu.fitPt(tau23mu::MuonFitType::TPFMS) > m_minMuPtCut)
            )
       )
      {
        event_.muons.push_back(ntupleMu);

        std::vector<Float_t> PTs = {ntupleMu.fitPt(tau23mu::MuonFitType::DEFAULT),
            ntupleMu.fitPt(tau23mu::MuonFitType::GLB),
            ntupleMu.fitPt(tau23mu::MuonFitType::TUNEP),
            ntupleMu.fitPt(tau23mu::MuonFitType::INNER),
            ntupleMu.fitPt(tau23mu::MuonFitType::PICKY),
            ntupleMu.fitPt(tau23mu::MuonFitType::DYT),
            ntupleMu.fitPt(tau23mu::MuonFitType::TPFMS)};
        eventId_.maxPTs.push_back(*std::max_element(PTs.begin(), PTs.end()));
      }

    }

    return event_.muons.size();

}

void Tau23MuNtupleMaker::fillTracks(const edm::Handle<std::vector<reco::Track> >& tracks,
      const edm::Handle<std::vector<reco::Vertex> >& vertexes,
      const edm::Handle<reco::BeamSpot>& beamSpot){

    std::vector<reco::Track>::const_iterator trIt  = tracks->begin();
    std::vector<reco::Track>::const_iterator trEnd = tracks->end();

    int nTrk = 0;
    for (; trIt != trEnd; ++trIt) 
    {

      const reco::Track tr = (*trIt);

      tau23mu::Track ntupleTr;

      ntupleTr.pt = tr.pt();
      ntupleTr.eta = tr.eta();
      ntupleTr.phi = tr.phi();
      ntupleTr.nChi2 = tr.normalizedChi2();
      ntupleTr.nValidHits = tr.numberOfValidHits();
      ntupleTr.charge = tr.charge();
      ntupleTr.dxy = tr.dxy();
      ntupleTr.dz = tr.dz();
      ntupleTr.vx = tr.vx();
      ntupleTr.vy = tr.vy();
      ntupleTr.vz = tr.vz();
      ntupleTr.dxyError = tr.dxyError();
      ntupleTr.dzError = tr.dzError();

      nTrk++;
      event_.tracks.push_back(ntupleTr);
    }
    event_.nTrk = nTrk;
}

bool Tau23MuNtupleMaker::getMuonChamberId(DetId & id, tau23mu::MuonDetType & det,
      Int_t & r, Int_t & phi, Int_t & eta) const
{

    if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT)
    {
      DTChamberId dtId(id.rawId());  

      det = tau23mu::MuonDetType::DT;
      r = dtId.station();
      phi = dtId.sector();
      eta = dtId.wheel();

      return true;
    }

    if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC)
    {
      CSCDetId cscId(id.rawId());

      det = tau23mu::MuonDetType::CSC;
      r = cscId.station() * cscId.zendcap();
      phi = cscId.chamber();
      eta = cscId.ring();

      return true;
    }

    return false;

}

void Tau23MuNtupleMaker::fillAnalysisTree(const edm::EventSetup & iSetup,
      const edm::Handle<edm::View<reco::Muon> > & muons,
      const edm::Handle<std::vector<reco::Track> >& tracks,
      const edm::Handle<std::vector<reco::Vertex> > & vertexes,
      const edm::Handle<std::vector<reco::Vertex> > & secVertexes,
      const edm::Handle<reco::BeamSpot> & beamSpot,
      const edm::Handle<edm::TriggerResults> & triggerResults, 
      const edm::Handle<trigger::TriggerEvent> & triggerEvent,
      const edm::TriggerNames & triggerNames,
      bool _runMC, bool do2mu_){

    tau23mu::triMuon_cand tmp_cand;
    tau23mu::diMuonTrk_cand tmp_2mutrk_cand;

    vector<size_t> goodMuonIndex(0);
    vector<size_t> muMatchedTrack(0);

    std::vector<reco::Track>::const_iterator trackIt = tracks->begin();
    std::vector<reco::Track>::const_iterator trackEnd = tracks->end();

    edm::View<reco::Muon>::const_iterator muonIt  = muons->begin();
    edm::View<reco::Muon>::const_iterator muonEnd = muons->end();

    Int_t mu_idx = -1;
    Int_t trk_idx = -1;
    Int_t n2muTrk = 0;
    Int_t n3mu = 0;

    for (; muonIt != muonEnd; ++muonIt) 
    {
      const reco::Muon& mu = (*muonIt);
      mu_idx++;
      if (mu.pt()<2 || abs(mu.eta())>2.4) continue;
      if (mu.isPFMuon() || mu.isGlobalMuon()) goodMuonIndex.push_back(mu_idx); // Store indices of all good muons
    }

    for (; trackIt != trackEnd; ++trackIt){
      const reco::Track& trk = (*trackIt);
      trk_idx++;

      for (; muonIt != muonEnd; ++muonIt){
        const reco::Muon& mu = (*muonIt);
        if (abs(trk.eta()-mu.eta())<0.001 && abs(trk.phi()-mu.phi())<0.001) muMatchedTrack.push_back(trk_idx);
      }
    }

    if (goodMuonIndex.size()<(do2mu_?2:3)) return; // return the number of good muons is less than 3(2) for Tau3Mu(DsPhiPi) analysis

    for (size_t i=0; i<goodMuonIndex.size()-1; ++i){
      const reco::Muon & mu1 = (*muons)[goodMuonIndex[i]]; // Select first good muon

      for (size_t j=i+1; j<goodMuonIndex.size(); ++j){
        const reco::Muon & mu2 = (*muons)[goodMuonIndex[j]]; // Select second good muon

        double dz_mu1mu2 = abs(mu1.innerTrack()->dz(beamSpot->position())-mu2.innerTrack()->dz(beamSpot->position()));
        if ( dz_mu1mu2 > dz_mumu_cut ) continue; // dz cut 

        double dR_mu1mu2 = deltaR(mu1.eta(), mu2.eta(), mu1.phi(), mu2.phi());
        if ( dR_mu1mu2 > dR_mumu_cut ) continue; // dR cut

        if ( goodMuonIndex.size() < j-1 ){

            for (size_t k=j+1; k<goodMuonIndex.size(); ++k){

              const reco::Muon & mu3 = (*muons)[goodMuonIndex[k]];
              size_t nPt2p5 = 0;
              if (mu1.pt()>2.5) nPt2p5++;
              if (mu2.pt()>2.5) nPt2p5++;
              if (mu3.pt()>2.5) nPt2p5++;

              if (nPt2p5<2) continue; // Require at least 2 muons with Pt more than 2.5 GeV

              double dz_mu2mu3 = abs(mu2.innerTrack()->dz(beamSpot->position())-mu3.innerTrack()->dz(beamSpot->position()));
              double dz_mu1mu3 = abs(mu1.innerTrack()->dz(beamSpot->position())-mu3.innerTrack()->dz(beamSpot->position()));

              if (dz_mu2mu3 > dz_mumu_cut || dz_mu1mu3 > dz_mumu_cut) continue; // dz cut

              double dR_mu2mu3 = deltaR(mu2.eta(), mu2.eta(), mu3.phi(), mu3.phi()); 
              double dR_mu1mu3 = deltaR(mu1.eta(), mu1.eta(), mu3.phi(), mu3.phi());

              if( dR_mu2mu3 > dR_mumu_cut || dR_mu1mu3 > dR_mumu_cut ) continue; // dR cut

              if ( mu1.charge()+mu2.charge()+mu3.charge() > 1) continue; // Require sum of charges to be +/- 1


              // Build tracks  and re-fit the vertex
              TrackRef trk1 = mu1.innerTrack();
              TrackRef trk2 = mu2.innerTrack();
              TrackRef trk3 = mu3.innerTrack();

              auto fv = fitTriMuVertex(iSetup, trk1, trk2, trk3);
              double fvnC2_tmp = fv.totalChiSquared()/fv.degreesOfFreedom();

              if (fvnC2_tmp > triMuFitVtx_nC2_cut ) continue; // Store all the good candidates

              tmp_cand.triMuFitVtx_nC2 = fvnC2_tmp;

              // Sort the muons by P
              vector<reco::Muon> sortedMu;

              sortedMu.push_back(tau23mu::maxPMu(tau23mu::maxPMu(mu1,mu2),mu3));
              sortedMu.push_back(tau23mu::minPMu(tau23mu::minPMu(maxPMu(mu1,mu2),tau23mu::maxPMu(mu2,mu3)),tau23mu::maxPMu(mu3,mu1)));
              sortedMu.push_back(tau23mu::minPMu(tau23mu::minPMu(mu1,mu2),mu3));

              TLorentzVector vtau;     
              TLorentzVector vec_mu1, vec_mu2, vec_mu3;

              vec_mu1.SetPtEtaPhiM(sortedMu[0].pt(), sortedMu[0].eta(), sortedMu[0].phi(), MU_MASS);
              vec_mu2.SetPtEtaPhiM(sortedMu[1].pt(), sortedMu[1].eta(), sortedMu[1].phi(), MU_MASS);
              vec_mu3.SetPtEtaPhiM(sortedMu[2].pt(), sortedMu[2].eta(), sortedMu[2].phi(), MU_MASS);

              vtau = vec_mu1+vec_mu2+vec_mu3;

              // Fill tree with trimuon information
              for (size_t it=0; it<3; ++it){

                tmp_cand.mu_p[it] = sortedMu[it].p();
                tmp_cand.mu_pt[it] = sortedMu[it].pt();
                tmp_cand.mu_pterr[it] = sortedMu[it].muonBestTrack()->ptError();
                tmp_cand.mu_bestPt[it] = sortedMu[it].muonBestTrack()->pt();
                tmp_cand.mu_eta[it] = sortedMu[it].eta();
                tmp_cand.mu_phi[it] = sortedMu[it].phi();
                tmp_cand.mu_charge[it] = sortedMu[it].charge();
                tmp_cand.mu_vx[it] = sortedMu[it].vx();
                tmp_cand.mu_vy[it] = sortedMu[it].vy();
                tmp_cand.mu_vz[it] = sortedMu[it].vz();

                // muon ID
                bool _isGlobal = sortedMu[it].isGlobalMuon();
                bool _isTracker = sortedMu[it].isTrackerMuon();
                bool _isTrackerArb = muon::isGoodMuon(sortedMu[it], muon::TrackerMuonArbitrated);
                bool _isRPC = sortedMu[it].isRPCMuon();
                bool _isStandAlone = sortedMu[it].isStandAloneMuon();
                bool _isPF = sortedMu[it].isPFMuon();

                bool _hasInnerTrack = !sortedMu[it].innerTrack().isNull();
                bool _hasTunePTrack = !sortedMu[it].tunePMuonBestTrack().isNull();
                bool _hasPickyTrack = !sortedMu[it].pickyTrack().isNull();
                bool _hasDytTrack = !sortedMu[it].dytTrack().isNull();
                bool _hasTpfmsTrack = !sortedMu[it].tpfmsTrack().isNull();

                // fill muon ID
                tmp_cand.mu_hasInnerTrack[it] = _hasInnerTrack;
                tmp_cand.mu_hasTunePTrack[it] = _hasTunePTrack;
                tmp_cand.mu_hasPickyTrack[it] = _hasPickyTrack;
                tmp_cand.mu_hasDytTrack[it] = _hasDytTrack;
                tmp_cand.mu_hasTpfmsTrack[it] = _hasTpfmsTrack;
                tmp_cand.mu_isGlobal[it] = _isGlobal;
                tmp_cand.mu_isTracker[it] = _isTracker;
                tmp_cand.mu_isTrackerArb[it] = _isTrackerArb;
                tmp_cand.mu_isRPC[it] = _isRPC;
                tmp_cand.mu_isStandAlone[it] = _isStandAlone;
                tmp_cand.mu_isPF[it] = _isPF;

                // Important variables
                tmp_cand.mu_uSta[it] = sortedMu[it].combinedQuality().updatedSta;
                tmp_cand.mu_trkKink[it] = sortedMu[it].combinedQuality().trkKink;
                tmp_cand.mu_log_glbKink[it] = TMath::Log(2+sortedMu[it].combinedQuality().glbKink);
                tmp_cand.mu_trkRelChi2[it] = sortedMu[it].combinedQuality().trkRelChi2;
                tmp_cand.mu_staRelChi2[it] = sortedMu[it].combinedQuality().staRelChi2;
                tmp_cand.mu_Chi2LP[it] = sortedMu[it].combinedQuality().chi2LocalPosition;
                tmp_cand.mu_Chi2LM[it] = sortedMu[it].combinedQuality().chi2LocalMomentum;
                tmp_cand.mu_localDist[it] = sortedMu[it].combinedQuality().localDistance;
                tmp_cand.mu_glbDEP[it] = sortedMu[it].combinedQuality().globalDeltaEtaPhi;
                tmp_cand.mu_tightMatch[it] = sortedMu[it].combinedQuality().tightMatch;
                tmp_cand.mu_glbTrkProb[it] = sortedMu[it].combinedQuality().glbTrackProbability;
                tmp_cand.mu_calEM[it] = sortedMu[it].calEnergy().em;
                tmp_cand.mu_calEMS9[it] = sortedMu[it].calEnergy().emS9;
                tmp_cand.mu_calEMS25[it] = sortedMu[it].calEnergy().emS25;
                tmp_cand.mu_calHad[it] = sortedMu[it].calEnergy().had;
                tmp_cand.mu_calHadS9[it] = sortedMu[it].calEnergy().hadS9;;
                tmp_cand.mu_nOMS[it] = sortedMu[it].numberOfMatchedStations();
                tmp_cand.mu_nOM[it] = sortedMu[it].numberOfMatches(reco::Muon::SegmentArbitration);
                tmp_cand.mu_comp2d[it] = muon::isGoodMuon(sortedMu[it], muon::TM2DCompatibilityTight);
                tmp_cand.mu_calocomp[it] = muon::caloCompatibility(sortedMu[it]);
                tmp_cand.mu_segcomp[it] = muon::segmentCompatibility(sortedMu[it]);

                // Muon outer track info
                tmp_cand.mu_oTrk_p[it] = _isStandAlone ? sortedMu[it].outerTrack()->p():-999;
                tmp_cand.mu_oTrk_pt[it] = _isStandAlone ? sortedMu[it].outerTrack()->pt():-999;
                tmp_cand.mu_oTrk_eta[it] = _isStandAlone ? sortedMu[it].outerTrack()->eta():-999;
                tmp_cand.mu_oTrk_phi[it] = _isStandAlone ? sortedMu[it].outerTrack()->phi(): -999;
                tmp_cand.mu_oTrk_nHits[it] = _isStandAlone ? sortedMu[it].outerTrack()->numberOfValidHits(): -999;

                tmp_cand.mu_oTrk_nC2[it] = _isStandAlone ? sortedMu[it].outerTrack()->normalizedChi2():-999;
                tmp_cand.mu_oTrk_MSWVH[it] = _isStandAlone ? sortedMu[it].outerTrack()->hitPattern().muonStationsWithValidHits():-999;
                tmp_cand.mu_oTrk_qprod[it] = _isStandAlone ? sortedMu[it].outerTrack()->charge()*sortedMu[it].innerTrack()->charge():-999;

                // Muon global track info
                tmp_cand.mu_glb_p[it] = _isGlobal ? sortedMu[it].globalTrack()->p():-999;
                tmp_cand.mu_glb_eta[it] = _isGlobal ? sortedMu[it].globalTrack()->eta():-999;
                tmp_cand.mu_glb_phi[it] = _isGlobal ? sortedMu[it].globalTrack()->phi():-999;
                tmp_cand.mu_glb_nC2[it] = _isGlobal ? sortedMu[it].globalTrack()->normalizedChi2():-999;
                tmp_cand.mu_glb_nOVMH[it] = _isGlobal ? sortedMu[it].globalTrack()->hitPattern().numberOfValidMuonHits():-999;
                tmp_cand.mu_glb_nHits[it] = _isGlobal ? sortedMu[it].globalTrack()->numberOfValidHits() : -999; 

                // muon inner track information
                tmp_cand.mu_inTrk_p[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->p():-999;
                tmp_cand.mu_inTrk_pt[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->pt():-999;
                tmp_cand.mu_inTrk_eta[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->eta():-999;
                tmp_cand.mu_inTrk_phi[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->phi():-999;
                tmp_cand.mu_inTrk_nC2[i] = _hasInnerTrack ? sortedMu[i].innerTrack()->normalizedChi2():-999;
                tmp_cand.mu_inTrk_trkLayWithMeas[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().trackerLayersWithMeasurement():-999;
                tmp_cand.mu_inTrk_pixLayWithMeas[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().pixelLayersWithMeasurement():-999;
                tmp_cand.mu_inTrk_nHitsTracker[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfValidTrackerHits():-999;
                tmp_cand.mu_inTrk_nHitsPixel[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfValidPixelHits():-999;
                tmp_cand.mu_inTrk_validFraction[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->validFraction():-999;
                tmp_cand.mu_inTrk_nLostTrkHits[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS):-999;
                tmp_cand.mu_inTrk_nLostTrkHits_in[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS):-999;
                tmp_cand.mu_inTrk_nLostTrkHits_out[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS):-999;
                tmp_cand.mu_inTrk_HP[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->quality(TrackBase::highPurity):-999;

                // Isolation variables
                reco::MuonIsolation IsoR03 = sortedMu[it].isolationR03();

                tmp_cand.mu_IsoR03_sumPt[it] = IsoR03.sumPt;
                tmp_cand.mu_IsoR03_nTrks[it] = IsoR03.nTracks;
                tmp_cand.mu_IsoR03_emEt[it] = IsoR03.emEt;
                tmp_cand.mu_IsoR03_hadEt[it] = IsoR03.hadEt;
                tmp_cand.mu_IsoR03_emVetoEt[it] = IsoR03.emVetoEt;
                tmp_cand.mu_IsoR03_hadVetoEt[it] = IsoR03.hadVetoEt;

                // PF isolation variables
                reco::MuonPFIsolation pfIso04 = sortedMu[it].pfIsolationR04();
                reco::MuonPFIsolation pfIso03 = sortedMu[it].pfIsolationR03();

                tmp_cand.chargedHadronIso[it] = pfIso04.sumChargedHadronPt;
                tmp_cand.chargedHadronIsoPU[it] = pfIso04.sumPUPt; 
                tmp_cand.neutralHadronIso[it] = pfIso04.sumNeutralHadronEt;
                tmp_cand.photonIso[it] = pfIso04.sumPhotonEt;

                tmp_cand.fits.push_back(tau23mu::MuonFit(sortedMu[it].pt(),sortedMu[it].eta(),sortedMu[it].phi(),
                        sortedMu[it].charge(),sortedMu[it].muonBestTrack()->ptError()));

                tmp_cand.fits.push_back(tau23mu::MuonFit(_hasInnerTrack ? sortedMu[it].innerTrack()->pt()  : -1000.,
                        _hasInnerTrack ? sortedMu[it].innerTrack()->eta() : -1000.,
                        _hasInnerTrack ? sortedMu[it].innerTrack()->phi() : -1000.,
                        _hasInnerTrack ? sortedMu[it].innerTrack()->charge()  : -1000.,
                        _hasInnerTrack ? sortedMu[it].innerTrack()->ptError() : -1000.));

                tmp_cand.fits.push_back(tau23mu::MuonFit(_isStandAlone ? sortedMu[it].outerTrack()->pt()  : -1000.,
                        _isStandAlone ? sortedMu[it].outerTrack()->eta() : -1000.,
                        _isStandAlone ? sortedMu[it].outerTrack()->phi() : -1000.,
                        _isStandAlone ? sortedMu[it].outerTrack()->charge()  : -1000.,
                        _isStandAlone ? sortedMu[it].outerTrack()->ptError() : -1000.));

                tmp_cand.fits.push_back(tau23mu::MuonFit(_isGlobal ? sortedMu[it].globalTrack()->pt()  : -1000.,
                        _isGlobal ? sortedMu[it].globalTrack()->eta() : -1000.,
                        _isGlobal ? sortedMu[it].globalTrack()->phi() : -1000.,
                        _isGlobal ? sortedMu[it].globalTrack()->charge()  : -1000.,
                        _isGlobal ? sortedMu[it].globalTrack()->ptError() : -1000.));

                tmp_cand.fits.push_back(tau23mu::MuonFit(_hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->pt()  : -1000.,
                        _hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->eta() : -1000.,
                        _hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->phi() : -1000.,
                        _hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->charge()  : -1000.,
                        _hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->ptError() : -1000.));

                tmp_cand.fits.push_back(tau23mu::MuonFit(_hasPickyTrack ? sortedMu[it].pickyTrack()->pt()  : -1000.,
                        _hasPickyTrack ? sortedMu[it].pickyTrack()->eta() : -1000.,
                        _hasPickyTrack ? sortedMu[it].pickyTrack()->phi() : -1000.,
                        _hasPickyTrack ? sortedMu[it].pickyTrack()->charge()  : -1000.,
                        _hasPickyTrack ? sortedMu[it].pickyTrack()->ptError() : -1000.));

                tmp_cand.fits.push_back(tau23mu::MuonFit(_hasDytTrack ? sortedMu[it].dytTrack()->pt()  : -1000.,
                        _hasDytTrack ? sortedMu[it].dytTrack()->eta() : -1000.,
                        _hasDytTrack ? sortedMu[it].dytTrack()->phi() : -1000.,
                        _hasDytTrack ? sortedMu[it].dytTrack()->charge()  : -1000.,
                        _hasDytTrack ? sortedMu[it].dytTrack()->ptError() : -1000.));

                tmp_cand.fits.push_back(tau23mu::MuonFit(_hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->pt()  : -1000.,
                        _hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->eta() : -1000.,
                        _hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->phi() : -1000.,
                        _hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->charge()  : -1000.,
                        _hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->ptError() : -1000.));

                tmp_cand.mu_dxyBest[it]  = -999; 
                tmp_cand.mu_dzBest[it] = -999; 

                tmp_cand.mu_dxybs[it] = _isGlobal ? sortedMu[it].globalTrack()->dxy(beamSpot->position()) :
                    _hasInnerTrack ? sortedMu[it].innerTrack()->dxy(beamSpot->position()) : -1000;
                tmp_cand.mu_dzbs[it]  = _isGlobal ? sortedMu[it].globalTrack()->dz(beamSpot->position()) :
                    _hasInnerTrack ? sortedMu[it].innerTrack()->dz(beamSpot->position()) : -1000;

                tmp_cand.mu_isSoft[it]  = 0; 
                tmp_cand.mu_isTight[it] = 0; 
                tmp_cand.mu_isHighPt[it]  = 0;
                tmp_cand.mu_isLoose[it] = muon::isLooseMuon(sortedMu[it])  ? 1 : 0; 
                tmp_cand.mu_isMedium[it]  = muon::isMediumMuon(sortedMu[it]) ? 1 : 0; 

                if ( sortedMu[it].isMatchesValid() && _isTrackerArb )
                {
                    for ( reco::MuonChamberMatch match : sortedMu[it].matches() )
                    {
                      tau23mu::ChambMatch ntupleMatch;

                      if ( getMuonChamberId(match.id,
                              ntupleMatch.type,ntupleMatch.r,
                              ntupleMatch.phi,ntupleMatch.eta)
                       )
                      {

                        ntupleMatch.x = sortedMu[it].trackX(match.station(),match.detector());
                        ntupleMatch.y = sortedMu[it].trackY(match.station(),match.detector());
                        ntupleMatch.dXdZ = sortedMu[it].trackDxDz(match.station(),match.detector());
                        ntupleMatch.dYdZ = sortedMu[it].trackDyDz(match.station(),match.detector());

                        ntupleMatch.errxTk = sortedMu[it].trackXErr(match.station(),match.detector());
                        ntupleMatch.erryTk = sortedMu[it].trackYErr(match.station(),match.detector());

                        ntupleMatch.errDxDzTk = sortedMu[it].trackDxDzErr(match.station(),match.detector());
                        ntupleMatch.errDyDzTk = sortedMu[it].trackDyDzErr(match.station(),match.detector());

                        ntupleMatch.dx = sortedMu[it].dX(match.station(),match.detector());
                        ntupleMatch.dy = sortedMu[it].dY(match.station(),match.detector());
                        ntupleMatch.dDxDz = sortedMu[it].dDxDz(match.station(),match.detector());
                        ntupleMatch.dDyDz = sortedMu[it].dDxDz(match.station(),match.detector());

                        ntupleMatch.errxSeg = sortedMu[it].segmentXErr(match.station(),match.detector());
                        ntupleMatch.errySeg = sortedMu[it].segmentYErr(match.station(),match.detector());
                        ntupleMatch.errDxDzSeg = sortedMu[it].segmentDxDzErr(match.station(),match.detector());
                        ntupleMatch.errDyDzSeg = sortedMu[it].segmentDyDzErr(match.station(),match.detector());

                        tmp_cand.matches.push_back(ntupleMatch);
                      }
                    }
                }
                tmp_cand.mu_isoPflow04[it] = (pfIso04.sumChargedHadronPt+ 
                      std::max(0.,pfIso04.sumPhotonEt+pfIso04.sumNeutralHadronEt - 0.5*pfIso04.sumPUPt)) / sortedMu[it].pt();

                tmp_cand.mu_isoPflow03[it] = (pfIso03.sumChargedHadronPt+ 
                      std::max(0.,pfIso03.sumPhotonEt+pfIso03.sumNeutralHadronEt - 0.5*pfIso03.sumPUPt)) / sortedMu[it].pt();

                if (vertexes->size() > 0) {
                    const reco::Vertex & vertex = vertexes->at(0);

                    tmp_cand.mu_dxy[it] = _isGlobal ? sortedMu[it].globalTrack()->dxy(vertex.position()) :
                      _hasInnerTrack ? sortedMu[it].innerTrack()->dxy(vertex.position()) : -1000;
                    tmp_cand.mu_dz[it] = _isGlobal ? sortedMu[it].globalTrack()->dz(vertex.position()) :
                      _hasInnerTrack ? sortedMu[it].innerTrack()->dz(vertex.position()) : -1000;

                    tmp_cand.mu_dxyBest[it]  = sortedMu[it].muonBestTrack()->dxy(vertex.position()); 
                    tmp_cand.mu_dzBest[it] = sortedMu[it].muonBestTrack()->dz(vertex.position()); 
                    tmp_cand.mu_dxyInner[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->dxy(vertex.position()):-999; 
                    tmp_cand.mu_dzInner[it]  = _hasInnerTrack ? sortedMu[it].innerTrack()->dz(vertex.position()):-999; 

                    tmp_cand.mu_isSoft[it]  = muon::isSoftMuon(sortedMu[it],vertex) ? 1 : 0; 
                    tmp_cand.mu_isTight[it] = muon::isTightMuon(sortedMu[it],vertex)  ? 1 : 0; 
                    tmp_cand.mu_isHighPt[it]  = muon::isHighPtMuon(sortedMu[it],vertex) ? 1 : 0;

                    tmp_cand.mu_dxy[it]  = _isGlobal ? sortedMu[it].globalTrack()->dxy(vertex.position()):-999;
                    tmp_cand.mu_dz[it]  = _isGlobal ? sortedMu[it].globalTrack()->dxy(vertex.position()):-999;

                }
                if(sortedMu[it].isTimeValid()) { 
                    tmp_cand.mu_TimeDof[it] = sortedMu[it].time().nDof; 
                    tmp_cand.mu_Time[it]  = sortedMu[it].time().timeAtIpInOut; 
                    tmp_cand.mu_TimeErr[it] = sortedMu[it].time().timeAtIpInOutErr; 
                } 
                else { 
                    tmp_cand.mu_TimeDof[it] = -999; 
                    tmp_cand.mu_Time[it]  = -999; 
                    tmp_cand.mu_TimeErr[it] = -999; 
                } 

                if(sortedMu[it].rpcTime().nDof > 0) { 
                    tmp_cand.mu_RpcTimeDof[it] = sortedMu[it].rpcTime().nDof; 
                    tmp_cand.mu_RpcTime[it]  = sortedMu[it].rpcTime().timeAtIpInOut; 
                    tmp_cand.mu_RpcTimeErr[it] = sortedMu[it].rpcTime().timeAtIpInOutErr; 
                } 
                else { 
                    tmp_cand.mu_RpcTimeDof[it] = -999; 
                    tmp_cand.mu_RpcTime[it] = -999; 
                    tmp_cand.mu_RpcTimeErr[it] = -999; 
                } 


                tmp_cand.mu_edxy[it] = _isGlobal ? sortedMu[it].globalTrack()->dxyError() : _hasInnerTrack ? sortedMu[it].innerTrack()->dxyError() : -1000;
                tmp_cand.mu_edz[it]  = _isGlobal ? sortedMu[it].globalTrack()->dzError()  : _hasInnerTrack ? sortedMu[it].innerTrack()->dzError() : -1000;

                tmp_cand.mu_TMOST[it] = muon::isGoodMuon(sortedMu[it], muon::TMOneStationTight);
                tmp_cand.mu_TMOSAT[it] = muon::isGoodMuon(sortedMu[it], muon::TMOneStationAngTight);
                tmp_cand.mu_TMLST[it] = muon::isGoodMuon(sortedMu[it], muon::TMLastStationTight);
                tmp_cand.mu_TMLSAT[it] = muon::isGoodMuon(sortedMu[it], muon::TMLastStationAngTight);
                tmp_cand.mu_TMLSOLPT[it] = muon::isGoodMuon(sortedMu[it], muon::TMLastStationOptimizedLowPtTight);
                tmp_cand.mu_TMLSOBLPT[it] = muon::isGoodMuon(sortedMu[it], muon::TMLastStationOptimizedBarrelLowPtTight);

                // Tau-Muon variables
                tmp_cand.mu_tau_dR[it] = deltaR(sortedMu[it].eta(), sortedMu[it].phi(), vtau.Eta(), vtau.Phi());

                // Good Global Muon, Tight Global Muon (previous definiitons)
                tmp_2mutrk_cand.mu_ggm[it] = (tmp_cand.mu_glb_nC2[it]<3 && tmp_cand.mu_Chi2LP[it]<12 && 
                      tmp_cand.mu_trkKink[it]<20 && tmp_cand.mu_segcomp[it]>0.303) ? 1:0;
                tmp_2mutrk_cand.mu_tgm[it] = (tmp_cand.mu_glb_nC2[it]<10 && tmp_cand.mu_isPF[it] && 
                      tmp_cand.mu_glb_nOVMH[it]>0 && tmp_cand.mu_nOMS[it]>1 && 
                      tmp_cand.mu_dxy[it]<0.2 && tmp_cand.mu_dz[it]<0.5 && 
                      tmp_cand.mu_nOVPH[it]>0 && tmp_cand.mu_inTrk_trkLayWithMeas[it]>5) ? 1:0; // d0 changed to dxy

              }

              tmp_cand.mu_pt_max =  TMath::Max(TMath::Max(mu1.pt(),mu2.pt()),mu3.pt()); //return highest pt muon 
              tmp_cand.mu_eta_max = TMath::Max(TMath::Max(mu1.eta(),mu2.eta()),mu3.eta()); //return highest eta muon 
              tmp_cand.mu_pt_max =  TMath::Min(TMath::Min(mu1.pt(),mu2.pt()),mu3.pt()); //return lowest pt muon 
              tmp_cand.mu_eta_max = TMath::Min(TMath::Min(mu1.eta(),mu2.eta()),mu3.eta()); //return lowest eta muon

              tmp_cand.mu_taudR_max = TMath::Max(TMath::Max(tmp_cand.mu_tau_dR[0],tmp_cand.mu_tau_dR[1]),tmp_cand.mu_tau_dR[2]);
              tmp_cand.mu_trkLayWithMeas_max = TMath::Max(TMath::Max(tmp_cand.mu_inTrk_trkLayWithMeas[0],tmp_cand.mu_inTrk_trkLayWithMeas[1]),tmp_cand.mu_inTrk_trkLayWithMeas[2]);

              // Re-assign and store dR and dz with sorted muons
              dR_mu1mu2 = deltaR(sortedMu[0].eta(), sortedMu[1].eta(), sortedMu[0].phi(), sortedMu[1].phi());
              dR_mu2mu3 = deltaR(sortedMu[1].eta(), sortedMu[2].eta(), sortedMu[1].phi(), sortedMu[2].phi());
              dR_mu1mu3 = deltaR(sortedMu[0].eta(), sortedMu[2].eta(), sortedMu[0].phi(), sortedMu[2].phi());

              tmp_cand.dR_mu1mu2 = dR_mu1mu2;
              tmp_cand.dR_mu2mu3 = dR_mu2mu3;
              tmp_cand.dR_mu1mu3 = dR_mu1mu3;

              tmp_cand.dz_mu1mu2 = abs(sortedMu[0].innerTrack()->dz(beamSpot->position())-sortedMu[1].innerTrack()->dz(beamSpot->position()));
              tmp_cand.dz_mu2mu3 = abs(sortedMu[1].innerTrack()->dz(beamSpot->position())-sortedMu[2].innerTrack()->dz(beamSpot->position()));
              tmp_cand.dz_mu1mu3 = abs(sortedMu[0].innerTrack()->dz(beamSpot->position())-sortedMu[2].innerTrack()->dz(beamSpot->position()));

              // Compute 3 mu invariant mass

              tmp_cand.M3mu = (vec_mu1+vec_mu2+vec_mu3).M();
              tmp_cand.M_mu1mu2 = (vec_mu1+vec_mu2).M();
              tmp_cand.M_mu2mu3 = (vec_mu2+vec_mu3).M();
              tmp_cand.M_mu1mu3 = (vec_mu1+vec_mu3).M();

              n3mu++;

              ////////////////
              // Refit tracks  
              ////////////////
              trk1 = sortedMu[0].innerTrack();
              trk2 = sortedMu[1].innerTrack();
              trk3 = sortedMu[2].innerTrack();
              ESHandle<TransientTrackBuilder> theB;
              iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
              vector<TransientTrack> t_trks;
              t_trks.push_back(theB->build(trk1));
              t_trks.push_back(theB->build(trk2));
              t_trks.push_back(theB->build(trk3)); 
              if(!fv.isValid()) { return; cout<<"[Tau23MuNtupleMaker]: Vertex Fit invalid!"<<endl; }

              // calculate closest approach ===> format has to change make a class and then get all the output results
              ClosestApproachInRPhi cApp12, cApp23, cApp31;
              cApp12.calculate(t_trks[0].initialFreeState(), t_trks[1].initialFreeState());
              cApp23.calculate(t_trks[1].initialFreeState(), t_trks[2].initialFreeState());
              cApp31.calculate(t_trks[2].initialFreeState(), t_trks[0].initialFreeState());

              if(!(cApp12.status()&&cApp23.status()&&cApp31.status())) { return; cout<<"[Tau23MuNtupleMaker]: DCA invalid!"<<endl; }
              tmp_cand.dca_mu1mu2 = cApp12.distance();
              tmp_cand.dca_mu2mu3 = cApp23.distance();
              tmp_cand.dca_mu1mu3 = cApp31.distance();
              tmp_cand.dca_max = TMath::Max(tmp_cand.dca_mu1mu2, TMath::Max(tmp_cand.dca_mu1mu3, tmp_cand.dca_mu2mu3));

              KalmanVertexFitter kvf(true);
              fv = kvf.vertex(t_trks);
              if(!fv.isValid()) { return; cout<<"[Tau23MuNtupleMaker]: Vertex Fit invalid!"<<endl; }

              TLorentzVector vtau_refit, vmu_refit;
              vtau_refit.SetPtEtaPhiM(0, 0, 0, 0);
              vector<TransientTrack>::const_iterator trkIt = fv.refittedTracks().begin();
              for(; trkIt != fv.refittedTracks().end(); ++ trkIt) {
                const reco::Track & trkrefit = trkIt->track();
                vmu_refit.SetPtEtaPhiM(trkrefit.pt(), trkrefit.eta(), trkrefit.phi(), MU_MASS);
                vtau_refit += vmu_refit;
              }

              tmp_cand.m3mu_refit = vtau_refit.M();

              tmp_cand.fv_tC2 = fv.totalChiSquared();
              tmp_cand.fv_dof = fv.degreesOfFreedom();
              tmp_cand.fv_nC2 = fv.totalChiSquared()/fv.degreesOfFreedom();
              tmp_cand.fv_Prob = TMath::Prob(tmp_cand.fv_tC2,(int)tmp_cand.fv_dof);

              vector<TransientTrack> t_trks12, t_trks23, t_trks31;
              t_trks12.push_back(theB->build(trk1)); t_trks12.push_back(theB->build(trk2));
              t_trks23.push_back(theB->build(trk2)); t_trks23.push_back(theB->build(trk3));
              t_trks31.push_back(theB->build(trk3)); t_trks31.push_back(theB->build(trk1));
              KalmanVertexFitter kvf_trks12, kvf_trks23, kvf_trks31;
              TransientVertex fv_trks12 = kvf_trks12.vertex(t_trks12);
              TransientVertex fv_trks23 = kvf_trks23.vertex(t_trks23);
              TransientVertex fv_trks31 = kvf_trks31.vertex(t_trks31);

              tmp_cand.fvwo_tC2[0] = fv_trks23.totalChiSquared();
              tmp_cand.fvwo_nC2[0] = tmp_cand.fvwo_tC2[0]/fv_trks23.degreesOfFreedom();
              tmp_cand.fvwo_tC2[1] = fv_trks31.totalChiSquared();
              tmp_cand.fvwo_nC2[1] = tmp_cand.fvwo_tC2[1]/fv_trks31.degreesOfFreedom();
              tmp_cand.fvwo_tC2[2] = fv_trks12.totalChiSquared();
              tmp_cand.fvwo_nC2[2] = tmp_cand.fvwo_tC2[2]/fv_trks12.degreesOfFreedom();

              TVector3 vtauxyz(vtau.Px(), vtau.Py(), vtau.Pz());

              ////////////////////
              // find the good PV
              ////////////////////
              Double_t PVZ = fv.position().z()-vtau.Pz()/vtau.Pt(); // original formula Double_t PVZ = fv.position().z()-fv_dxy*vtau.Pz()/vtau.Pt(); with no fv_dxy defined ???
              Double_t dispv1 = 999, dispvgen=999;
              tmp_cand.dphi_pv = -999;
              //int ipvPVZ = 99, ipvgen = 99;
              size_t ipv1, ipv2, ipv_gen;
              double gen_pv = 0;

              for(size_t jpv = 0; jpv < vertexes->size(); jpv++) {
                const Vertex & vi = (*vertexes)[jpv];

                if(abs(vi.position().Z()-PVZ)<dispv1){
                    dispv1=abs(vi.position().Z()-PVZ);
                    ipv1=jpv;
                }
                if(abs(vi.position().Z()-gen_pv)<dispvgen){
                    dispvgen=abs(vi.position().Z()-gen_pv);
                    ipv_gen=jpv;
                }

                TVector3 dv_3d(fv.position().x() - vi.x(), fv.position().y() - vi.y(), fv.position().z() - vi.z());
                Double_t cos_dphi_3d = dv_3d.Dot(vtauxyz)/(dv_3d.Mag()*vtauxyz.Mag());
                if(cos_dphi_3d > tmp_cand.dphi_pv){
                    tmp_cand.dphi_pv = cos_dphi_3d;
                    ipv2=jpv;
                }

              }
              const Vertex & pv0 = (*vertexes)[ipv2];

              //////////////////////////////////////////////////
              // refit PV with and w.o. the 3 mu
              //////////////////////////////////////////////////

              tmp_cand.pv1_tC2 = 999; 
              tmp_cand.pv1_nC2 = 999; 
              tmp_cand.pv2_tC2 = 999; 
              tmp_cand.pv2_nC2 = 999;

              vector<TransientTrack> pv_trks;
              TransientVertex pv2, pv1;

              iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
              for(Vertex::trackRef_iterator itk = pv0.tracks_begin(); itk != pv0.tracks_end(); itk++) {
                if((**itk).pt()>1) {
                    if(deltaR(sortedMu[0].eta(), sortedMu[0].phi(), (**itk).eta(), (**itk).phi())<0.01)continue;
                    if(deltaR(sortedMu[1].eta(), sortedMu[1].phi(), (**itk).eta(), (**itk).phi())<0.01)continue;
                    if(deltaR(sortedMu[2].eta(), sortedMu[2].phi(), (**itk).eta(), (**itk).phi())<0.01)continue;
                }
                pv_trks.push_back(theB->build(**itk));
              }

              if(pv_trks.size()>1) {
                KalmanVertexFitter kvf_pv;
                pv1 = kvf_pv.vertex(pv_trks);
                if(pv1.isValid()){
                    tmp_cand.pv1_tC2 = pv1.totalChiSquared();
                    tmp_cand.pv1_nC2 = pv1.totalChiSquared()/pv1.degreesOfFreedom();
                }

                // adding the 3 mu tracks
                pv_trks.push_back(theB->build(trk1));
                pv_trks.push_back(theB->build(trk2));
                pv_trks.push_back(theB->build(trk3));
                pv2 = kvf_pv.vertex(pv_trks);
                if(pv2.isValid()){
                    tmp_cand.pv2_tC2 = pv2.totalChiSquared();
                    tmp_cand.pv2_nC2 = pv2.totalChiSquared()/pv2.degreesOfFreedom();
                }
              }

              Vertex pvv = pv0;  // the final PV
              if(pv1.isValid()) pvv = Vertex(pv1);
              math::XYZPoint pv1P = math::XYZPoint(pvv.x(), pvv.y(), pvv.z());

              tmp_cand.d0[0] = abs(sortedMu[0].innerTrack()->dxy(pv1P));
              tmp_cand.d0[1] = abs(sortedMu[1].innerTrack()->dxy(pv1P));
              tmp_cand.d0[2] = abs(sortedMu[2].innerTrack()->dxy(pv1P));

              for (size_t _i=0; _i<3; ++_i) tmp_cand.d0sig[_i]=-999;
              GlobalVector dir1(sortedMu[0].px(),sortedMu[0].py(), sortedMu[0].pz());
              GlobalVector dir2(sortedMu[1].px(), sortedMu[1].py(), sortedMu[1].pz());
              GlobalVector dir3(sortedMu[2].px(), sortedMu[2].py(), sortedMu[2].pz());
              std::pair<bool, Measurement1D> ip2d_1 = IPTools::signedTransverseImpactParameter(t_trks[0], dir1, pvv);
              std::pair<bool, Measurement1D> ip2d_2 = IPTools::signedTransverseImpactParameter(t_trks[1], dir2, pvv);
              std::pair<bool, Measurement1D> ip2d_3 = IPTools::signedTransverseImpactParameter(t_trks[2], dir3, pvv);
              if(ip2d_1.first) tmp_cand.d0sig[0] = abs(ip2d_1.second.value()/ip2d_1.second.error());
              if(ip2d_2.first) tmp_cand.d0sig[1] = abs(ip2d_2.second.value()/ip2d_2.second.error());
              if(ip2d_3.first) tmp_cand.d0sig[2] = abs(ip2d_3.second.value()/ip2d_3.second.error());

              ////////////////////
              // displacement 2D
              TVector3 dv_2d(fv.position().x() - pv1P.x(), fv.position().y() - pv1P.y(), 0);
              TVector3 vtauxy(vtau.Px(), vtau.Py(), 0);
              tmp_cand.fv_cosdphi = dv_2d.Dot(vtauxy)/(dv_2d.Perp()*vtauxy.Perp());
              VertexDistanceXY vdistXY;
              Measurement1D distXY = vdistXY.distance(Vertex(fv), pvv);
              tmp_cand.fv_dxy = distXY.value();
              tmp_cand.fv_dxysig = distXY.significance();
              tmp_cand.fv_ppdl3D = distXY.value() * tmp_cand.fv_cosdphi * tmp_cand.M3mu/vtauxy.Perp();

              ////////////////////
              // displacement 3D
              TVector3 dv_3d(fv.position().x() - pv1P.x(), fv.position().y() - pv1P.y(), fv.position().z() - pv1P.z());
              //TVector3 vtauxyz(vtau.Px(), vtau.Py(), vtau.Pz());
              tmp_cand.fv_cosdphi3D = dv_3d.Dot(vtauxyz)/(dv_3d.Mag()*vtauxyz.Mag());
              VertexDistance3D dist;
              tmp_cand.fv_d3D = dist.distance(Vertex(fv), pvv).value(); // = dv_reco.Mag() ??
              tmp_cand.fv_d3Dsig = dist.distance(Vertex(fv), pvv).significance();
              tmp_cand.fv_ppdl3D = tmp_cand.fv_d3D*tmp_cand.fv_cosdphi3D*tmp_cand.M3mu/vtau.P();

              vector<double> softmueta, softmuphi;

              tmp_cand.pv_nSoftMu = 0; // count number of muons from the primary vertex

              ////////////////////
              // Soft Muon Info
              ////////////////////
              for(size_t _softMu = 0; _softMu < muons->size(); ++ _softMu) {
                if(_softMu==goodMuonIndex[i])continue;
                if(_softMu==goodMuonIndex[j])continue;
                if(_softMu==goodMuonIndex[k])continue;
                const reco::Muon & m_1 = (*muons)[i];
                if(!(abs(m_1.eta())<2.4)) continue;
                if(!(muon::isGoodMuon(m_1, muon::TMOneStationTight))) continue;
                if(!(m_1.innerTrack()->hitPattern().trackerLayersWithMeasurement()>5))continue;
                if(!(m_1.innerTrack()->hitPattern().pixelLayersWithMeasurement()>0))continue;
                //if(!(abs(m_1.innerTrack()->dxy(pv0.position())) < .3))continue;
                if(!(abs(m_1.innerTrack()->dz(pv1P)) < 1))continue;
                tmp_cand.pv_nSoftMu++;
                softmueta.push_back(m_1.eta());
                softmuphi.push_back(m_1.phi());
              }

              ////////////////////
              // secondary vertices
              ////////////////////
              tmp_cand.n_sv = 0;

              for(size_t isv = 0; isv < secVertexes->size(); isv++) {
                const Vertex & sv = (*secVertexes)[isv];
                if(abs(sv.p4().M()-KPM_MASS)<K_mass_cut && sv.tracksSize()==2) continue; // no Ks

                double dx = sv.x()-pv1P.x();
                double dy = sv.y()-pv1P.y();
                double dz = sv.z()-pv1P.z();

                TVector3 sv_reco(dx, dy, dz);
                TVector3 svxyz(sv.p4().Px(), sv.p4().Py(), sv.p4().Pz());

                VertexDistance3D distsv;

                auto temp_sv_d3D = distsv.distance(sv, pvv).value();
                auto temp_sv_cosdphi_3D = sv_reco.Dot(svxyz)/(sv_reco.Mag()*svxyz.Mag());

                tmp_cand.sv_cosdphi3D.push_back(temp_sv_cosdphi_3D);
                tmp_cand.sv_d3D.push_back(temp_sv_d3D);
                tmp_cand.sv_overlap.push_back(deltaR(sv_reco.Eta(), sv_reco.Phi(), dv_3d.Eta(), dv_3d.Phi()));
                tmp_cand.sv_d3Dsig.push_back(distsv.distance(sv, pvv).significance());
                tmp_cand.sv_ppdl3D.push_back(temp_sv_d3D*temp_sv_cosdphi_3D*sv.p4().M()/sv.p4().P());

                Int_t temp_sv_nmu =0;

                for(Vertex::trackRef_iterator itk = sv.tracks_begin(); itk != sv.tracks_end(); itk++) {
                    for(size_t imu = 0; imu < softmueta.size(); imu ++) {
                      if(deltaR(softmueta[imu], softmuphi[imu], (**itk).eta(), (**itk).phi())<0.01)
                        temp_sv_nmu ++;
                    }
                }

                tmp_cand.sv_nmu.push_back(temp_sv_nmu);
                tmp_cand.sv_mass.push_back(sv.p4().M());
                tmp_cand.sv_pt.push_back(sv.p4().Pt());
                tmp_cand.sv_dz.push_back(abs(dz));
                tmp_cand.sv_ntrk.push_back(sv.tracksSize());
                tmp_cand.n_sv++;
              }

              tmp_cand.dzpv[0] = 1;
              tmp_cand.dzpv[1] = 1;
              tmp_cand.dzpv[2] = 1;

              tmp_cand.dzpv[0] = abs(sortedMu[0].innerTrack()->dz(pv1P));
              tmp_cand.dzpv[1] = abs(sortedMu[1].innerTrack()->dz(pv1P));
              tmp_cand.dzpv[2] = abs(sortedMu[2].innerTrack()->dz(pv1P));

              for(size_t jpv = 0; jpv < vertexes->size(); jpv++) {
                if(jpv==ipv2)continue;
                const Vertex & vi = (*vertexes)[jpv];
                if(abs(sortedMu[0].innerTrack()->dz(vi.position()))<tmp_cand.mu_dz[0]) tmp_cand.dzpv[0]=-1;
                if(abs(sortedMu[1].innerTrack()->dz(vi.position()))<tmp_cand.mu_dz[1]) tmp_cand.dzpv[1]=-1;
                if(abs(sortedMu[2].innerTrack()->dz(vi.position()))<tmp_cand.mu_dz[2]) tmp_cand.dzpv[2]=-1;
              }

              ////////////////////
              // Track Isolation
              // How to decide if a track is associated with a certain PV ?
              ////////////////////

              double pttrk_tau = 0, pttrk_tau0p5 = 0,  pttrk_m1 = 0, pttrk_m2 = 0, pttrk_m3 = 0;

              tmp_cand.mindca_iso = 999;
              tmp_cand.mindca_iso05 = 999;

              // Initialize track parameters
              tmp_cand.ntrk_tau = 0; 
              tmp_cand.ntrk_tau0p5 = 0; 
              tmp_cand.ntrk_tau_b = 0; 
              tmp_cand.ntrk_sum = 0;
              for (size_t it = 0; it<3; ++it) tmp_cand.ntrk[it] = 0;  
              tmp_cand.ntrk0p1 = 0; 
              tmp_cand.ntrk0p2 = 0;
              tmp_cand.ntrk0p5 = 0;
              tmp_cand.maxdxy_pv0 = -1;

              math::XYZPoint fvP = math::XYZPoint(fv.position().x(), fv.position().y(), fv.position().z());

              for(size_t itk = 0; itk < tracks->size(); itk++) {
                const reco::Track & t = (*tracks)[itk];
                if(!(t.quality(TrackBase::tight)))continue;
                if(deltaR(sortedMu[0].eta(), sortedMu[0].phi(), t.eta(), t.phi())<0.01)continue;
                if(deltaR(sortedMu[1].eta(), sortedMu[1].phi(), t.eta(), t.phi())<0.01)continue;
                if(deltaR(sortedMu[2].eta(), sortedMu[2].phi(), t.eta(), t.phi())<0.01)continue;

                double dz = abs(t.dz(fvP));
                double dxy = abs(t.dxy(fvP));
                double dca_fv = sqrt(dz*dz+dxy*dxy);
                double dr_tau = deltaR(t.eta(), t.phi(), vtau.Eta(), vtau.Phi());

                // iso no. 1b - using pt_min, drtau_max of the 3 mu
                if(t.pt() > 0.33*tmp_cand.mu_pt_min && dr_tau < 3.*tmp_cand.mu_taudR_max && dca_fv<0.05 ) {
                    pttrk_tau += t.pt();
                    tmp_cand.ntrk_tau++; // iso 3b
                    if(dca_fv<tmp_cand.mindca_iso) tmp_cand.mindca_iso=dca_fv; // iso 4b
                } 

                if(t.pt()<1.0) continue;  // was 1.2
                // iso no. 1
                if(dr_tau < 0.5 && dca_fv<0.05 ) {
                    pttrk_tau0p5 += t.pt();
                    tmp_cand.ntrk_tau0p5++; // iso 3
                    //if(dca_fv<mindca_iso05)mindca_iso05=dca_fv; // iso 4
                }

                if(dca_fv<0.05) tmp_cand.ntrk_tau_b++; // iso 3b
                if(dca_fv<tmp_cand.mindca_iso05) tmp_cand.mindca_iso05=dca_fv; // iso 4

                TransientTrack trkiso = theB->build(t);
                ClosestApproachInRPhi cAppm1, cAppm2, cAppm3;
                cAppm1.calculate(trkiso.initialFreeState(), t_trks[0].initialFreeState());
                cAppm2.calculate(trkiso.initialFreeState(), t_trks[1].initialFreeState());
                cAppm3.calculate(trkiso.initialFreeState(), t_trks[2].initialFreeState());
                if(!(cAppm1.status()&&cAppm2.status()&&cAppm3.status())) continue;

                // iso no. 2
                if(deltaR(t.eta(), t.phi(), sortedMu[0].eta(), sortedMu[0].phi()) < 0.3 && cAppm1.distance() < 0.1) {// && dz1 < .3) 
                    tmp_cand.ntrk[0]++;
                    pttrk_m1 += t.pt();
                }
                if(deltaR(t.eta(), t.phi(), sortedMu[1].eta(), sortedMu[1].phi()) < 0.3 && cAppm2.distance() < 0.1) {//&& dz2 < .3) 
                    tmp_cand.ntrk[1]++;
                    pttrk_m2 += t.pt();
                }
                if(deltaR(t.eta(), t.phi(), sortedMu[2].eta(), sortedMu[2].phi()) < 0.3 && cAppm3.distance() < 0.1) {//&& dz3 < .3) 
                    tmp_cand.ntrk[2]++;
                    pttrk_m3 += t.pt();
                }
                if( (deltaR(t.eta(), t.phi(), sortedMu[0].eta(), sortedMu[0].phi()) < 0.3 && cAppm1.distance() < 0.1 )
                      ||(deltaR(t.eta(), t.phi(), sortedMu[1].eta(), sortedMu[1].phi()) < 0.3 && cAppm2.distance() < 0.1 )
                      ||(deltaR(t.eta(), t.phi(), sortedMu[2].eta(), sortedMu[2].phi()) < 0.3 && cAppm3.distance() < 0.1 )
                  ) tmp_cand.ntrk_sum++;


                // displaced track counting
                // only tracks consistent with PV
                double dz_pv0=abs(t.dz(pv1P));
                if(!(dz_pv0 < 1))continue;
                double dxy_pv0 = abs(t.dxy(pv1P));
                if(dxy_pv0>0.1) tmp_cand.ntrk0p1++;
                if(dxy_pv0>0.2) tmp_cand.ntrk0p2++;
                if(dxy_pv0>0.5) tmp_cand.ntrk0p5++;
                if(dxy_pv0>tmp_cand.maxdxy_pv0) tmp_cand.maxdxy_pv0 = dxy_pv0;

                tmp_cand.trkrel_tau = pttrk_tau/vtau.Pt();
                tmp_cand.trkrel_tau0p5 = pttrk_tau0p5/vtau.Pt();
                tmp_cand.trkrel[0] = pttrk_m1/sortedMu[0].pt(); 
                tmp_cand.trkrel[1] = pttrk_m2/sortedMu[1].pt();
                tmp_cand.trkrel[2] = pttrk_m3/sortedMu[2].pt();
                tmp_cand.trkrel_max = TMath::Max(tmp_cand.trkrel[0], TMath::Max(tmp_cand.trkrel[1], tmp_cand.trkrel[2]));

              }
              event_.triMuon_cand_coll.push_back(tmp_cand);
              n3mu++;
            }

            std::vector<HitInfo> hits; // where to put hit??

            //***************************************************
            // Find 2mu+1track candidates
            //***************************************************

            if (do2mu_ && mu1.pt()>2.5 && mu2.pt()>2.5){
				  double min_fvnC_2mu1tk = 10;
              size_t trk_idx=-1, tmp_iter=0;
              int j1 =-1, j2=-1;

              std::vector<reco::Track>::const_iterator trkIt = tracks->begin();
              std::vector<reco::Track>::const_iterator trkEnd = tracks->end();

              // sort muons in the order of momentum
              reco::Muon sortedMu[2];
              sortedMu[0] = mu1.p()>mu2.p() ? mu1:mu2; 
              sortedMu[1] = mu1.p()<mu2.p() ? mu1:mu2;

              for (; trkIt!=trkEnd; ++trkIt){

                trk_idx++;

                if (trk_idx==muMatchedTrack[tmp_iter]) tmp_iter++; // skip if the track is matched to a muon

                else {
                    const reco::Track& trk = (*trkIt);

                    if (abs(mu1.charge()+mu2.charge()+trk.charge())>1) continue;

                    // Build tracks  and re-fit the vertex
                    TrackRef trk1 = mu1.innerTrack();
                    TrackRef trk2 = mu2.innerTrack();
                    TrackRef trk3 = TrackRef(tracks, trk_idx);

                    TransientVertex fv = fitTriMuVertex(iSetup, trk1, trk2, trk3);
                    if(!fv.isValid()) continue;
                    double fv_tC2 = fv.totalChiSquared();
                    double fv_dOF = fv.degreesOfFreedom();
                    double fv_nC = fv_tC2/fv_dOF;
                    if(fv_nC>5) continue;
                    if(fv_nC < min_fvnC_2mu1tk){

                      // iTrk=itk; don't know what this stands for

                      if(sortedMu[0].p()>sortedMu[1].p()){
                        j1=goodMuonIndex[i]; j2=goodMuonIndex[j];
                      }
                      else {j1=goodMuonIndex[j]; j2=goodMuonIndex[i];}
                      min_fvnC_2mu1tk = fv_nC;
                    }

                    if(!fv.isValid()) continue;

                    double fvnC2_tmp = fv.totalChiSquared()/fv.degreesOfFreedom();

                    if (fvnC2_tmp > diMuTrkFitVtx_nC2_cut ) continue; // Store all the good candidates

                    TLorentzVector vec_mu1, vec_mu2, vec_trk, vds;

                    vec_mu1.SetPtEtaPhiM(sortedMu[0].pt(), sortedMu[0].eta(), sortedMu[0].phi(), MU_MASS);
                    vec_mu2.SetPtEtaPhiM(sortedMu[1].pt(), sortedMu[1].eta(), sortedMu[1].phi(), MU_MASS);
                    vec_trk.SetPtEtaPhiM(trk.pt(), trk.eta(), trk.phi(), PIPM_MASS);

                    vds = vec_mu1+vec_mu2+vec_trk;

                    // Fill 2 muon info
                    for (size_t it=0; it<2; ++it){

                      tmp_2mutrk_cand.mu_p[it] = sortedMu[it].p();
                      tmp_2mutrk_cand.mu_pt[it] = sortedMu[it].pt();
                      tmp_2mutrk_cand.mu_pterr[it] = sortedMu[it].muonBestTrack()->ptError();
                      tmp_2mutrk_cand.mu_bestPt[it] = sortedMu[it].muonBestTrack()->pt();
                      tmp_2mutrk_cand.mu_eta[it] = sortedMu[it].eta();
                      tmp_2mutrk_cand.mu_phi[it] = sortedMu[it].phi();
                      tmp_2mutrk_cand.mu_charge[it] = sortedMu[it].charge();
                      tmp_2mutrk_cand.mu_vx[it] = sortedMu[it].vx();
                      tmp_2mutrk_cand.mu_vy[it] = sortedMu[it].vy();
                      tmp_2mutrk_cand.mu_vz[it] = sortedMu[it].vz();

                      // muon ID
                      bool _isGlobal = sortedMu[it].isGlobalMuon();
                      bool _isTracker = sortedMu[it].isTrackerMuon();
                      bool _isTrackerArb = muon::isGoodMuon(sortedMu[it], muon::TrackerMuonArbitrated);
                      bool _isRPC = sortedMu[it].isRPCMuon();
                      bool _isStandAlone = sortedMu[it].isStandAloneMuon();
                      bool _isPF = sortedMu[it].isPFMuon();

                      bool _hasInnerTrack = !sortedMu[it].innerTrack().isNull();
                      bool _hasTunePTrack = !sortedMu[it].tunePMuonBestTrack().isNull();
                      bool _hasPickyTrack = !sortedMu[it].pickyTrack().isNull();
                      bool _hasDytTrack = !sortedMu[it].dytTrack().isNull();
                      bool _hasTpfmsTrack = !sortedMu[it].tpfmsTrack().isNull();

                      // fill muon ID
                      tmp_2mutrk_cand.mu_hasInnerTrack[it] = _hasInnerTrack;
                      tmp_2mutrk_cand.mu_hasTunePTrack[it] = _hasTunePTrack;
                      tmp_2mutrk_cand.mu_hasPickyTrack[it] = _hasPickyTrack;
                      tmp_2mutrk_cand.mu_hasDytTrack[it] = _hasDytTrack;
                      tmp_2mutrk_cand.mu_hasTpfmsTrack[it] = _hasTpfmsTrack;
                      tmp_2mutrk_cand.mu_isGlobal[it] = _isGlobal;
                      tmp_2mutrk_cand.mu_isTracker[it] = _isTracker;
                      tmp_2mutrk_cand.mu_isTrackerArb[it] = _isTrackerArb;
                      tmp_2mutrk_cand.mu_isRPC[it] = _isRPC;
                      tmp_2mutrk_cand.mu_isStandAlone[it] = _isStandAlone;
                      tmp_2mutrk_cand.mu_isPF[it] = _isPF; 

                      tmp_2mutrk_cand.mu_uSta[it] = sortedMu[it].combinedQuality().updatedSta;
                      tmp_2mutrk_cand.mu_trkKink[it] = sortedMu[it].combinedQuality().trkKink;
                      tmp_2mutrk_cand.mu_log_glbKink[it] = TMath::Log(2+sortedMu[it].combinedQuality().glbKink);
                      tmp_2mutrk_cand.mu_trkRelChi2[it] = sortedMu[it].combinedQuality().trkRelChi2;
                      tmp_2mutrk_cand.mu_staRelChi2[it] = sortedMu[it].combinedQuality().staRelChi2;
                      tmp_2mutrk_cand.mu_Chi2LP[it] = sortedMu[it].combinedQuality().chi2LocalPosition;
                      tmp_2mutrk_cand.mu_Chi2LM[it] = sortedMu[it].combinedQuality().chi2LocalMomentum;
                      tmp_2mutrk_cand.mu_localDist[it] = sortedMu[it].combinedQuality().localDistance;
                      tmp_2mutrk_cand.mu_glbDEP[it] = sortedMu[it].combinedQuality().globalDeltaEtaPhi;
                      tmp_2mutrk_cand.mu_tightMatch[it] = sortedMu[it].combinedQuality().tightMatch;
                      tmp_2mutrk_cand.mu_glbTrkProb[it] = sortedMu[it].combinedQuality().glbTrackProbability;
                      tmp_2mutrk_cand.mu_calEM[it] = sortedMu[it].calEnergy().em;
                      tmp_2mutrk_cand.mu_calEMS9[it] = sortedMu[it].calEnergy().emS9;
                      tmp_2mutrk_cand.mu_calEMS25[it] = sortedMu[it].calEnergy().emS25;
                      tmp_2mutrk_cand.mu_calHad[it] = sortedMu[it].calEnergy().had;
                      tmp_2mutrk_cand.mu_calHadS9[it] = sortedMu[it].calEnergy().hadS9;;
                      tmp_2mutrk_cand.mu_nOMS[it] = sortedMu[it].numberOfMatchedStations();
                      tmp_2mutrk_cand.mu_nOM[it] = sortedMu[it].numberOfMatches(reco::Muon::SegmentArbitration);
                      tmp_2mutrk_cand.mu_comp2d[it] = muon::isGoodMuon(sortedMu[it], muon::TM2DCompatibilityTight);
                      tmp_2mutrk_cand.mu_calocomp[it] = muon::caloCompatibility(sortedMu[it]);
                      tmp_2mutrk_cand.mu_segcomp[it] = muon::segmentCompatibility(sortedMu[it]);

                      // Muon outer track info
                      tmp_2mutrk_cand.mu_oTrk_p[it] = _isStandAlone ? sortedMu[it].outerTrack()->p():-999;
                      tmp_2mutrk_cand.mu_oTrk_eta[it] = _isStandAlone ? sortedMu[it].outerTrack()->eta():-999;
                      tmp_2mutrk_cand.mu_oTrk_phi[it] = _isStandAlone ? sortedMu[it].outerTrack()->phi(): -999;
                      tmp_2mutrk_cand.mu_oTrk_nC2[it] = _isStandAlone ? sortedMu[it].outerTrack()->normalizedChi2():-999;
                      tmp_2mutrk_cand.mu_oTrk_MSWVH[it] = _isStandAlone ? sortedMu[it].outerTrack()->hitPattern().muonStationsWithValidHits():-999;
                      tmp_2mutrk_cand.mu_oTrk_qprod[it] = _isStandAlone ? sortedMu[it].outerTrack()->charge()*sortedMu[it].innerTrack()->charge():-999;

                      // Muon global track info
                      tmp_2mutrk_cand.mu_glb_p[it] = _isGlobal ? sortedMu[it].globalTrack()->p():-999;
                      tmp_2mutrk_cand.mu_glb_eta[it] = _isGlobal ? sortedMu[it].globalTrack()->eta():-999;
                      tmp_2mutrk_cand.mu_glb_phi[it] = _isGlobal ? sortedMu[it].globalTrack()->phi():-999;
                      tmp_2mutrk_cand.mu_glb_nC2[it] = _isGlobal ? sortedMu[it].globalTrack()->normalizedChi2():-999;
                      tmp_2mutrk_cand.mu_glb_nOVMH[it] = _isGlobal ? sortedMu[it].globalTrack()->hitPattern().numberOfValidMuonHits():-999;

                      // muon inner track information
                      tmp_2mutrk_cand.mu_inTrk_p[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->p():-999;
                      tmp_2mutrk_cand.mu_inTrk_pt[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->pt():-999;
                      tmp_2mutrk_cand.mu_inTrk_eta[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->eta():-999;
                      tmp_2mutrk_cand.mu_inTrk_phi[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->phi():-999;
                      tmp_2mutrk_cand.mu_inTrk_nC2[i] = _hasInnerTrack ? sortedMu[i].innerTrack()->normalizedChi2():-999;
                      tmp_2mutrk_cand.mu_inTrk_trkLayWithMeas[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().trackerLayersWithMeasurement():-999;
                      tmp_2mutrk_cand.mu_inTrk_pixLayWithMeas[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().pixelLayersWithMeasurement():-999;
                      tmp_2mutrk_cand.mu_inTrk_nHitsTracker[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfValidTrackerHits():-999;
                      tmp_2mutrk_cand.mu_inTrk_nHitsPixel[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfValidTrackerHits():-999;
                      tmp_2mutrk_cand.mu_inTrk_validFraction[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->validFraction():-999;
                      tmp_2mutrk_cand.mu_inTrk_nLostTrkHits[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS):-999;
                      tmp_2mutrk_cand.mu_inTrk_nLostTrkHits_in[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS):-999;
                      tmp_2mutrk_cand.mu_inTrk_nLostTrkHits_out[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS):-999;
                      tmp_2mutrk_cand.mu_inTrk_HP[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->quality(TrackBase::highPurity):-999;

                      // Isolation variables
                      reco::MuonIsolation IsoR03 = sortedMu[it].isolationR03();

                      tmp_2mutrk_cand.mu_IsoR03_sumPt[it] = IsoR03.sumPt;
                      tmp_2mutrk_cand.mu_IsoR03_nTrks[it] = IsoR03.nTracks;
                      tmp_2mutrk_cand.mu_IsoR03_emEt[it] = IsoR03.emEt;
                      tmp_2mutrk_cand.mu_IsoR03_hadEt[it] = IsoR03.hadEt;
                      tmp_2mutrk_cand.mu_IsoR03_emVetoEt[it] = IsoR03.emVetoEt;
                      tmp_2mutrk_cand.mu_IsoR03_hadVetoEt[it] = IsoR03.hadVetoEt;

                      // PF isolation variables
                      reco::MuonPFIsolation pfIso04 = sortedMu[it].pfIsolationR04();
                      reco::MuonPFIsolation pfIso03 = sortedMu[it].pfIsolationR03();

                      tmp_2mutrk_cand.chargedHadronIso[it] = pfIso04.sumChargedHadronPt;
                      tmp_2mutrk_cand.chargedHadronIsoPU[it] = pfIso04.sumPUPt; 
                      tmp_2mutrk_cand.neutralHadronIso[it] = pfIso04.sumNeutralHadronEt;
                      tmp_2mutrk_cand.photonIso[it] = pfIso04.sumPhotonEt;

                      tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(sortedMu[it].pt(),sortedMu[it].eta(),sortedMu[it].phi(),
                              sortedMu[it].charge(),sortedMu[it].muonBestTrack()->ptError()));

                      tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(_hasInnerTrack ? sortedMu[it].innerTrack()->pt()  : -1000.,
                              _hasInnerTrack ? sortedMu[it].innerTrack()->eta() : -1000.,
                              _hasInnerTrack ? sortedMu[it].innerTrack()->phi() : -1000.,
                              _hasInnerTrack ? sortedMu[it].innerTrack()->charge()  : -1000.,
                              _hasInnerTrack ? sortedMu[it].innerTrack()->ptError() : -1000.));

                      tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(_isStandAlone ? sortedMu[it].outerTrack()->pt()  : -1000.,
                              _isStandAlone ? sortedMu[it].outerTrack()->eta() : -1000.,
                              _isStandAlone ? sortedMu[it].outerTrack()->phi() : -1000.,
                              _isStandAlone ? sortedMu[it].outerTrack()->charge()  : -1000.,
                              _isStandAlone ? sortedMu[it].outerTrack()->ptError() : -1000.));

                      tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(_isGlobal ? sortedMu[it].globalTrack()->pt()  : -1000.,
                              _isGlobal ? sortedMu[it].globalTrack()->eta() : -1000.,
                              _isGlobal ? sortedMu[it].globalTrack()->phi() : -1000.,
                              _isGlobal ? sortedMu[it].globalTrack()->charge()  : -1000.,
                              _isGlobal ? sortedMu[it].globalTrack()->ptError() : -1000.));

                      tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(_hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->pt()  : -1000.,
                              _hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->eta() : -1000.,
                              _hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->phi() : -1000.,
                              _hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->charge()  : -1000.,
                              _hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->ptError() : -1000.));

                      tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(_hasPickyTrack ? sortedMu[it].pickyTrack()->pt()  : -1000.,
                              _hasPickyTrack ? sortedMu[it].pickyTrack()->eta() : -1000.,
                              _hasPickyTrack ? sortedMu[it].pickyTrack()->phi() : -1000.,
                              _hasPickyTrack ? sortedMu[it].pickyTrack()->charge()  : -1000.,
                              _hasPickyTrack ? sortedMu[it].pickyTrack()->ptError() : -1000.));

                      tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(_hasDytTrack ? sortedMu[it].dytTrack()->pt()  : -1000.,
                              _hasDytTrack ? sortedMu[it].dytTrack()->eta() : -1000.,
                              _hasDytTrack ? sortedMu[it].dytTrack()->phi() : -1000.,
                              _hasDytTrack ? sortedMu[it].dytTrack()->charge()  : -1000.,
                              _hasDytTrack ? sortedMu[it].dytTrack()->ptError() : -1000.));

                      tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(_hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->pt()  : -1000.,
                              _hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->eta() : -1000.,
                              _hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->phi() : -1000.,
                              _hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->charge()  : -1000.,
                              _hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->ptError() : -1000.));

                      tmp_2mutrk_cand.mu_dxyBest[it]  = -999; 
                      tmp_2mutrk_cand.mu_dzBest[it] = -999; 
                      tmp_2mutrk_cand.mu_dxyInner[it] = -999; 
                      tmp_2mutrk_cand.mu_dzInner[it]  = -999; 

                      tmp_2mutrk_cand.mu_dxybs[it] = _isGlobal ? sortedMu[it].globalTrack()->dxy(beamSpot->position()) :
                        _hasInnerTrack ? sortedMu[it].innerTrack()->dxy(beamSpot->position()) : -1000;
                      tmp_2mutrk_cand.mu_dzbs[it]  = _isGlobal ? sortedMu[it].globalTrack()->dz(beamSpot->position()) :
                        _hasInnerTrack ? sortedMu[it].innerTrack()->dz(beamSpot->position()) : -1000;

                      tmp_2mutrk_cand.mu_isSoft[it]  = 0; 
                      tmp_2mutrk_cand.mu_isTight[it] = 0; 
                      tmp_2mutrk_cand.mu_isHighPt[it]  = 0;
                      tmp_2mutrk_cand.mu_isLoose[it] = muon::isLooseMuon(sortedMu[it])  ? 1 : 0; 
                      tmp_2mutrk_cand.mu_isMedium[it]  = muon::isMediumMuon(sortedMu[it]) ? 1 : 0; 
                      if ( sortedMu[it].isMatchesValid() && _isTrackerArb )
                      {
                        for ( reco::MuonChamberMatch match : sortedMu[it].matches() )
                        {
                            tau23mu::ChambMatch ntupleMatch;

                            if ( getMuonChamberId(match.id,
                                    ntupleMatch.type,ntupleMatch.r,
                                    ntupleMatch.phi,ntupleMatch.eta)
                             )
                            {

                              ntupleMatch.x = sortedMu[it].trackX(match.station(),match.detector());
                              ntupleMatch.y = sortedMu[it].trackY(match.station(),match.detector());
                              ntupleMatch.dXdZ = sortedMu[it].trackDxDz(match.station(),match.detector());
                              ntupleMatch.dYdZ = sortedMu[it].trackDyDz(match.station(),match.detector());

                              ntupleMatch.errxTk = sortedMu[it].trackXErr(match.station(),match.detector());
                              ntupleMatch.erryTk = sortedMu[it].trackYErr(match.station(),match.detector());

                              ntupleMatch.errDxDzTk = sortedMu[it].trackDxDzErr(match.station(),match.detector());
                              ntupleMatch.errDyDzTk = sortedMu[it].trackDyDzErr(match.station(),match.detector());

                              ntupleMatch.dx = sortedMu[it].dX(match.station(),match.detector());
                              ntupleMatch.dy = sortedMu[it].dY(match.station(),match.detector());
                              ntupleMatch.dDxDz = sortedMu[it].dDxDz(match.station(),match.detector());
                              ntupleMatch.dDyDz = sortedMu[it].dDxDz(match.station(),match.detector());

                              ntupleMatch.errxSeg = sortedMu[it].segmentXErr(match.station(),match.detector());
                              ntupleMatch.errySeg = sortedMu[it].segmentYErr(match.station(),match.detector());
                              ntupleMatch.errDxDzSeg = sortedMu[it].segmentDxDzErr(match.station(),match.detector());
                              ntupleMatch.errDyDzSeg = sortedMu[it].segmentDyDzErr(match.station(),match.detector());

                              tmp_2mutrk_cand.matches.push_back(ntupleMatch);
                            }
                        }
                      }
                      tmp_2mutrk_cand.mu_isoPflow04[it] = (pfIso04.sumChargedHadronPt+ 
                            std::max(0.,pfIso04.sumPhotonEt+pfIso04.sumNeutralHadronEt - 0.5*pfIso04.sumPUPt)) / sortedMu[it].pt();

                      tmp_2mutrk_cand.mu_isoPflow03[it] = (pfIso03.sumChargedHadronPt+ 
                            std::max(0.,pfIso03.sumPhotonEt+pfIso03.sumNeutralHadronEt - 0.5*pfIso03.sumPUPt)) / sortedMu[it].pt();

                      if (vertexes->size() > 0) {
                        const reco::Vertex & vertex = vertexes->at(0);

                        tmp_2mutrk_cand.mu_dxy[it] = _isGlobal ? sortedMu[it].globalTrack()->dxy(vertex.position()) :
                            _hasInnerTrack ? sortedMu[it].innerTrack()->dxy(vertex.position()) : -1000;
                        tmp_2mutrk_cand.mu_dz[it] = _isGlobal ? sortedMu[it].globalTrack()->dz(vertex.position()) :
                            _hasInnerTrack ? sortedMu[it].innerTrack()->dz(vertex.position()) : -1000;

                        tmp_2mutrk_cand.mu_dxyBest[it]  = sortedMu[it].muonBestTrack()->dxy(vertex.position()); 
                        tmp_2mutrk_cand.mu_dzBest[it] = sortedMu[it].muonBestTrack()->dz(vertex.position()); 
                        if(_hasInnerTrack) { 
                            tmp_2mutrk_cand.mu_dxyInner[it] = sortedMu[it].innerTrack()->dxy(vertex.position()); 
                            tmp_2mutrk_cand.mu_dzInner[it]  = sortedMu[it].innerTrack()->dz(vertex.position()); 
                        } 

                        tmp_2mutrk_cand.mu_isSoft[it]  = muon::isSoftMuon(sortedMu[it],vertex) ? 1 : 0; 
                        tmp_2mutrk_cand.mu_isTight[it] = muon::isTightMuon(sortedMu[it],vertex)  ? 1 : 0; 
                        tmp_2mutrk_cand.mu_isHighPt[it]  = muon::isHighPtMuon(sortedMu[it],vertex) ? 1 : 0;

                        tmp_2mutrk_cand.mu_dxy[it]  = _isGlobal ? sortedMu[it].globalTrack()->dxy(vertex.position()):-999;
                        tmp_2mutrk_cand.mu_dz[it]  = _isGlobal ? sortedMu[it].globalTrack()->dxy(vertex.position()):-999;

                      }
                      if(sortedMu[it].isTimeValid()) { 
                        tmp_2mutrk_cand.mu_TimeDof[it] = sortedMu[it].time().nDof; 
                        tmp_2mutrk_cand.mu_Time[it]  = sortedMu[it].time().timeAtIpInOut; 
                        tmp_2mutrk_cand.mu_TimeErr[it] = sortedMu[it].time().timeAtIpInOutErr; 
                      } 
                      else { 
                        tmp_2mutrk_cand.mu_TimeDof[it] = -999; 
                        tmp_2mutrk_cand.mu_Time[it]  = -999; 
                        tmp_2mutrk_cand.mu_TimeErr[it] = -999; 
                      } 

                      if(sortedMu[it].rpcTime().nDof > 0) { 
                        tmp_2mutrk_cand.mu_RpcTimeDof[it] = sortedMu[it].rpcTime().nDof; 
                        tmp_2mutrk_cand.mu_RpcTime[it]  = sortedMu[it].rpcTime().timeAtIpInOut; 
                        tmp_2mutrk_cand.mu_RpcTimeErr[it] = sortedMu[it].rpcTime().timeAtIpInOutErr; 
                      } 
                      else { 
                        tmp_2mutrk_cand.mu_RpcTimeDof[it] = -999; 
                        tmp_2mutrk_cand.mu_RpcTime[it] = -999; 
                        tmp_2mutrk_cand.mu_RpcTimeErr[it] = -999; 
                      } 



                      tmp_2mutrk_cand.mu_edxy[it] = _isGlobal ? sortedMu[it].globalTrack()->dxyError() : _hasInnerTrack ? sortedMu[it].innerTrack()->dxyError() : -1000;
                      tmp_2mutrk_cand.mu_edz[it]  = _isGlobal ? sortedMu[it].globalTrack()->dzError()  : _hasInnerTrack ? sortedMu[it].innerTrack()->dzError() : -1000;

                      tmp_2mutrk_cand.mu_TMOST[it] = muon::isGoodMuon(sortedMu[it], muon::TMOneStationTight);
                      tmp_2mutrk_cand.mu_TMOSAT[it] = muon::isGoodMuon(sortedMu[it], muon::TMOneStationAngTight);
                      tmp_2mutrk_cand.mu_TMLST[it] = muon::isGoodMuon(sortedMu[it], muon::TMLastStationTight);
                      tmp_2mutrk_cand.mu_TMLSAT[it] = muon::isGoodMuon(sortedMu[it], muon::TMLastStationAngTight);
                      tmp_2mutrk_cand.mu_TMLSOLPT[it] = muon::isGoodMuon(sortedMu[it], muon::TMLastStationOptimizedLowPtTight);
                      tmp_2mutrk_cand.mu_TMLSOBLPT[it] = muon::isGoodMuon(sortedMu[it], muon::TMLastStationOptimizedBarrelLowPtTight);

                      // Tau-Muon variables
                      tmp_2mutrk_cand.mu_ds_dR[it] = deltaR(sortedMu[it].eta(), sortedMu[it].phi(), vds.Eta(), vds.Phi());

                      // Good Global Muon, Tight Global Muon (previous definiitons)
                      tmp_2mutrk_cand.mu_ggm[it] = (tmp_2mutrk_cand.mu_glb_nC2[it]<3 && tmp_2mutrk_cand.mu_Chi2LP[it]<12 && 
                            tmp_2mutrk_cand.mu_trkKink[it]<20 && tmp_2mutrk_cand.mu_segcomp[it]>0.303) ? 1:0;
                      tmp_2mutrk_cand.mu_tgm[it] = (tmp_2mutrk_cand.mu_glb_nC2[it]<10 && tmp_2mutrk_cand.mu_isPF[it] && 
                            tmp_2mutrk_cand.mu_glb_nOVMH[it]>0 && tmp_2mutrk_cand.mu_nOMS[it]>1 && 
                            tmp_2mutrk_cand.mu_dxy[it]<0.2 && tmp_2mutrk_cand.mu_dz[it]<0.5 && 
                            tmp_2mutrk_cand.mu_nOVPH[it]>0 && tmp_2mutrk_cand.mu_inTrk_trkLayWithMeas[it]>5) ? 1:0; // d0 changed to dxy
                    }

                    tmp_2mutrk_cand.mu_dsdR_max = TMath::Max(TMath::Max(tmp_2mutrk_cand.mu_ds_dR[0],tmp_2mutrk_cand.mu_ds_dR[1]),tmp_2mutrk_cand.mu_ds_dR[2]);
                    tmp_2mutrk_cand.mu_trkLayWithMeas_max = TMath::Max(TMath::Max(tmp_2mutrk_cand.mu_inTrk_trkLayWithMeas[0],tmp_2mutrk_cand.mu_inTrk_trkLayWithMeas[1]),tmp_2mutrk_cand.mu_inTrk_trkLayWithMeas[2]);

                    // Fill track info
                    tmp_2mutrk_cand.trk_pterr = trk.ptError();
                    tmp_2mutrk_cand.trk_trkhp = trk.quality(TrackBase::highPurity);
                    tmp_2mutrk_cand.trk_charge = trk.charge();
                    tmp_2mutrk_cand.trk_p = trk.p();
                    tmp_2mutrk_cand.trk_pt = trk.pt();
                    tmp_2mutrk_cand.trk_eta = trk.eta();
                    tmp_2mutrk_cand.trk_phi = trk.phi();
                    tmp_2mutrk_cand.trk_nOVPH = trk.hitPattern().numberOfValidPixelHits();
                    tmp_2mutrk_cand.trk_iTvF = trk.validFraction();
                    tmp_2mutrk_cand.trk_tLWM = trk.hitPattern().trackerLayersWithMeasurement();
                    tmp_2mutrk_cand.trk_pLWM = trk.hitPattern().pixelLayersWithMeasurement();
                    tmp_2mutrk_cand.trk_nOVTH = trk.hitPattern().numberOfValidTrackerHits();
                    tmp_2mutrk_cand.trk_nOLTH = trk.hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS);
                    tmp_2mutrk_cand.trk_nOLTHin = trk.hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);
                    tmp_2mutrk_cand.trk_nOLTHout = trk.hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS);
                    tmp_2mutrk_cand.trk_iTnC = trk.normalizedChi2();

                    // Compute 2mu + trk invaiant mass
                    tmp_2mutrk_cand.M2muTrk = (vec_mu1+vec_mu2+vec_trk).M();
                    tmp_2mutrk_cand.diMuTrkFitVtx_nC2 = fvnC2_tmp;

                    tmp_2mutrk_cand.M_mu1mu2 = (vec_mu1+vec_mu2).M();
                    tmp_2mutrk_cand.M_mu1trk = (vec_mu1+vec_trk).M();
                    tmp_2mutrk_cand.M_mu2trk = (vec_mu2+vec_trk).M();


                    // Refit tracks  
                    trk1 = sortedMu[0].innerTrack();
                    trk2 = sortedMu[1].innerTrack();
                    ESHandle<TransientTrackBuilder> theB;
                    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
                    vector<TransientTrack> t_trks;

                    t_trks.push_back(theB->build(trk1));
                    t_trks.push_back(theB->build(trk2));
                    t_trks.push_back(theB->build(trk)); 

                    if(!fv.isValid()) { return; cout<<"[Tau23MuNtupleMaker]: Vertex Fit invalid!"<<endl; }

                    // calculate closest approach ===> format has to change make a class and then get all the output results
                    ClosestApproachInRPhi cApp12, cApp23, cApp31;
                    cApp12.calculate(t_trks[0].initialFreeState(), t_trks[1].initialFreeState());
                    cApp23.calculate(t_trks[1].initialFreeState(), t_trks[2].initialFreeState());
                    cApp31.calculate(t_trks[2].initialFreeState(), t_trks[0].initialFreeState());

                    if(!(cApp12.status()&&cApp23.status()&&cApp31.status())) { return; cout<<"[Tau23MuNtupleMaker]: DCA invalid!"<<endl; }
                    tmp_2mutrk_cand.dca_mu1mu2 = cApp12.distance();
                    tmp_2mutrk_cand.dca_mu2trk = cApp23.distance();
                    tmp_2mutrk_cand.dca_mu1trk = cApp31.distance();
                    tmp_2mutrk_cand.dca_max = TMath::Max(tmp_2mutrk_cand.dca_mu1mu2, TMath::Max(tmp_2mutrk_cand.dca_mu1trk, tmp_2mutrk_cand.dca_mu2trk));

                    KalmanVertexFitter kvf(true);
                    fv = kvf.vertex(t_trks);
                    if(!fv.isValid()) { return; cout<<"[Tau23MuNtupleMaker]: Vertex Fit unvalid!"<<endl; }

                    TLorentzVector vds_refit, vmu_refit;
                    vds_refit.SetPtEtaPhiM(0, 0, 0, 0);
                    vector<TransientTrack>::const_iterator trkIt = fv.refittedTracks().begin();
                    for(; trkIt != fv.refittedTracks().end(); ++ trkIt) {
                      const reco::Track & trkrefit = trkIt->track();
                      vmu_refit.SetPtEtaPhiM(trkrefit.pt(), trkrefit.eta(), trkrefit.phi(), PIPM_MASS);
                      vds_refit += vmu_refit;
                    }

                    tmp_2mutrk_cand.m2muTrk_refit = vds_refit.M();

                    tmp_2mutrk_cand.fv_tC2 = fv.totalChiSquared();
                    tmp_2mutrk_cand.fv_dof = fv.degreesOfFreedom();
                    tmp_2mutrk_cand.fv_nC2 = fv.totalChiSquared()/fv.degreesOfFreedom();
                    tmp_2mutrk_cand.fv_Prob = TMath::Prob(tmp_2mutrk_cand.fv_tC2,(int)tmp_2mutrk_cand.fv_dof);

                    vector<TransientTrack> t_trks12, t_trks23, t_trks31;
                    t_trks12.push_back(theB->build(trk1)); t_trks12.push_back(theB->build(trk2));
                    t_trks23.push_back(theB->build(trk2)); t_trks23.push_back(theB->build(trk));
                    t_trks31.push_back(theB->build(trk)); t_trks31.push_back(theB->build(trk1));
                    KalmanVertexFitter kvf_trks12, kvf_trks23, kvf_trks31;
                    TransientVertex fv_trks12 = kvf_trks12.vertex(t_trks12);
                    TransientVertex fv_trks23 = kvf_trks23.vertex(t_trks23);
                    TransientVertex fv_trks31 = kvf_trks31.vertex(t_trks31);

                    tmp_2mutrk_cand.fvwo_tC2[0] = fv_trks23.totalChiSquared();
                    tmp_2mutrk_cand.fvwo_nC2[0] = tmp_2mutrk_cand.fvwo_tC2[0]/fv_trks23.degreesOfFreedom();
                    tmp_2mutrk_cand.fvwo_tC2[1] = fv_trks31.totalChiSquared();
                    tmp_2mutrk_cand.fvwo_nC2[1] = tmp_2mutrk_cand.fvwo_tC2[1]/fv_trks31.degreesOfFreedom();
                    tmp_2mutrk_cand.fvwo_tC2[2] = fv_trks12.totalChiSquared();
                    tmp_2mutrk_cand.fvwo_nC2[2] = tmp_2mutrk_cand.fvwo_tC2[2]/fv_trks12.degreesOfFreedom();

                    TVector3 vdsxyz(vds.Px(), vds.Py(), vds.Pz());
                    ////////////////////
                    // find the good PV
                    ////////////////////
                    Double_t PVZ = fv.position().z()-vds.Pz()/vds.Pt();
                    Double_t dispv1 = 999, dispvgen=999;
                    tmp_2mutrk_cand.dphi_pv = -999;
                    //int ipvPVZ = 99, ipvgen = 99;
                    size_t ipv1, ipv2, ipv_gen;
                    double gen_pv = 0; // need to change the gen_pv for MC 

                    for(size_t jpv = 0; jpv < vertexes->size(); jpv++) {
                      const Vertex & vi = (*vertexes)[jpv];

                      if(abs(vi.position().Z()-PVZ)<dispv1){
                        dispv1=abs(vi.position().Z()-PVZ);
                        ipv1=jpv;
                      }
                      if(abs(vi.position().Z()-gen_pv)<dispvgen){
                        dispvgen=abs(vi.position().Z()-gen_pv);
                        ipv_gen=jpv;
                      }

                      TVector3 dv_3d(fv.position().x() - vi.x(), fv.position().y() - vi.y(), fv.position().z() - vi.z());
                      Double_t cos_dphi_3d = dv_3d.Dot(vdsxyz)/(dv_3d.Mag()*vdsxyz.Mag());
                      if(cos_dphi_3d > tmp_2mutrk_cand.dphi_pv){
                        tmp_2mutrk_cand.dphi_pv = cos_dphi_3d;
                        ipv2=jpv;
                      }

                    }
                    const Vertex & pv0 = (*vertexes)[ipv2];

                    //////////////////////////////////////////////////
                    // refit PV with and w.o. the 3 mu
                    //////////////////////////////////////////////////

                    tmp_2mutrk_cand.pv1_tC2 = 999; 
                    tmp_2mutrk_cand.pv1_nC2 = 999; 
                    tmp_2mutrk_cand.pv2_tC2 = 999; 
                    tmp_2mutrk_cand.pv2_nC2 = 999;

                    vector<TransientTrack> pv_trks;
                    TransientVertex pv2, pv1;

                    // Fix the track in both classes 
                    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
                    for(Vertex::trackRef_iterator itk = pv0.tracks_begin(); itk != pv0.tracks_end(); itk++) {
                      if((**itk).pt()>1) {
                        if (deltaR(sortedMu[0].eta(), sortedMu[1].phi(), (**itk).eta(), (**itk).phi())<0.01) continue;
                        if (deltaR(sortedMu[1].eta(), sortedMu[1].phi(), (**itk).eta(), (**itk).phi())<0.01) continue;
                        if(deltaR((**itk).eta(), (**itk).phi(), trk.eta(), trk.phi())<0.01)continue;
                      }
                      pv_trks.push_back(theB->build(**itk));
                    }

                    if(pv_trks.size()>1) {
                      KalmanVertexFitter kvf_pv;
                      pv1 = kvf_pv.vertex(pv_trks);
                      if(pv1.isValid()){
                        tmp_2mutrk_cand.pv1_tC2 = pv1.totalChiSquared();
                        tmp_2mutrk_cand.pv1_nC2 = pv1.totalChiSquared()/pv1.degreesOfFreedom();
                      }

                      // adding the 3 mu tracks
                      pv_trks.push_back(theB->build(trk1));
                      pv_trks.push_back(theB->build(trk2));
                      pv_trks.push_back(theB->build(trk3));
                      pv2 = kvf_pv.vertex(pv_trks);
                      if(pv2.isValid()){
                        tmp_2mutrk_cand.pv2_tC2 = pv2.totalChiSquared();
                        tmp_2mutrk_cand.pv2_nC2 = pv2.totalChiSquared()/pv2.degreesOfFreedom();
                      }
                    }

                    Vertex pvv = pv0;  // the final PV
                    if(pv1.isValid()) pvv = Vertex(pv1);
                    math::XYZPoint pv1P = math::XYZPoint(pvv.x(), pvv.y(), pvv.z());

                    tmp_2mutrk_cand.d0[0] = abs(sortedMu[0].innerTrack()->dxy(pv1P));
                    tmp_2mutrk_cand.d0[1] = abs(sortedMu[1].innerTrack()->dxy(pv1P));
                    tmp_2mutrk_cand.d0[2] = abs(trk.dxy(pv1P));

                    for (size_t _i=0; _i<3; ++_i) tmp_2mutrk_cand.d0sig[_i]=-999;
                    GlobalVector dir1(sortedMu[0].px(),sortedMu[0].py(), sortedMu[0].pz());
                    GlobalVector dir2(sortedMu[1].px(), sortedMu[1].py(), sortedMu[1].pz());
                    GlobalVector dir3(trk.px(), trk.py(), trk.pz());
                    std::pair<bool, Measurement1D> ip2d_1 = IPTools::signedTransverseImpactParameter(t_trks[0], dir1, pvv);
                    std::pair<bool, Measurement1D> ip2d_2 = IPTools::signedTransverseImpactParameter(t_trks[1], dir2, pvv);
                    std::pair<bool, Measurement1D> ip2d_3 = IPTools::signedTransverseImpactParameter(t_trks[2], dir3, pvv);
                    if(ip2d_1.first) tmp_2mutrk_cand.d0sig[0] = abs(ip2d_1.second.value()/ip2d_1.second.error());
                    if(ip2d_2.first) tmp_2mutrk_cand.d0sig[1] = abs(ip2d_2.second.value()/ip2d_2.second.error());
                    if(ip2d_3.first) tmp_2mutrk_cand.d0sig[2] = abs(ip2d_3.second.value()/ip2d_3.second.error());

                    ////////////////////
                    // displacement 2D
                    TVector3 dv_2d(fv.position().x() - pv1P.x(), fv.position().y() - pv1P.y(), 0);
                    TVector3 vdsxy(vds.Px(), vds.Py(), 0);
                    tmp_2mutrk_cand.fv_cosdphi = dv_2d.Dot(vdsxy)/(dv_2d.Perp()*vdsxy.Perp());
                    VertexDistanceXY vdistXY;
                    Measurement1D distXY = vdistXY.distance(Vertex(fv), pvv);
                    tmp_2mutrk_cand.fv_dxy = distXY.value();
                    tmp_2mutrk_cand.fv_dxysig = distXY.significance();
                    tmp_2mutrk_cand.fv_ppdl3D = distXY.value() * tmp_2mutrk_cand.fv_cosdphi * tmp_2mutrk_cand.M2muTrk/vdsxy.Perp();

                    ////////////////////
                    // displacement 3D
                    TVector3 dv_3d(fv.position().x() - pv1P.x(), fv.position().y() - pv1P.y(), fv.position().z() - pv1P.z());
                    //TVector3 vdsxyz(vds.Px(), vds.Py(), vds.Pz());
                    tmp_2mutrk_cand.fv_cosdphi3D = dv_3d.Dot(vdsxyz)/(dv_3d.Mag()*vdsxyz.Mag());
                    VertexDistance3D dist;
                    tmp_2mutrk_cand.fv_d3D = dist.distance(Vertex(fv), pvv).value(); // = dv_reco.Mag() ??
                    tmp_2mutrk_cand.fv_d3Dsig = dist.distance(Vertex(fv), pvv).significance();
                    tmp_2mutrk_cand.fv_ppdl3D = tmp_2mutrk_cand.fv_d3D*tmp_2mutrk_cand.fv_cosdphi3D*tmp_2mutrk_cand.M2muTrk/vds.P();

                    vector<double> softmueta, softmuphi;

                    tmp_2mutrk_cand.pv_nSoftMu = 0; // count number of muons from the primary vertex

                    ////////////////////
                    // Soft Muon Info
                    ////////////////////
                    for(size_t _softMu = 0; _softMu < muons->size(); ++ _softMu) {
                      if(_softMu==goodMuonIndex[i])continue;
                      if(_softMu==goodMuonIndex[j])continue;
                      const reco::Muon & m_1 = (*muons)[i];
                      if(!(abs(m_1.eta())<2.4)) continue;
                      if(!(muon::isGoodMuon(m_1, muon::TMOneStationTight))) continue;
                      if(!(m_1.innerTrack()->hitPattern().trackerLayersWithMeasurement()>5))continue;
                      if(!(m_1.innerTrack()->hitPattern().pixelLayersWithMeasurement()>0))continue;
                      //if(!(abs(m_1.innerTrack()->dxy(pv0.position())) < .3))continue;
                      if(!(abs(m_1.innerTrack()->dz(pv1P)) < 1))continue;
                      tmp_2mutrk_cand.pv_nSoftMu++;
                      softmueta.push_back(m_1.eta());
                      softmuphi.push_back(m_1.phi());
                    }

                    ////////////////////
                    // secondary vertices
                    ////////////////////
                    tmp_2mutrk_cand.n_sv = 0;

                    for(size_t isv = 0; isv < secVertexes->size(); isv++) {
                      const Vertex & sv = (*secVertexes)[isv];
                      if(abs(sv.p4().M()-KPM_MASS)<K_mass_cut && sv.tracksSize()==2) continue; // no Ks

                      double dx = sv.x()-pv1P.x();
                      double dy = sv.y()-pv1P.y();
                      double dz = sv.z()-pv1P.z();

                      TVector3 sv_reco(dx, dy, dz);
                      TVector3 svxyz(sv.p4().Px(), sv.p4().Py(), sv.p4().Pz());

                      VertexDistance3D distsv;

                      auto temp_sv_d3D = distsv.distance(sv, pvv).value();
                      auto temp_sv_cosdphi_3D = sv_reco.Dot(svxyz)/(sv_reco.Mag()*svxyz.Mag());

                      tmp_2mutrk_cand.sv_cosdphi3D.push_back(temp_sv_cosdphi_3D);
                      tmp_2mutrk_cand.sv_d3D.push_back(temp_sv_d3D);
                      tmp_2mutrk_cand.sv_overlap.push_back(deltaR(sv_reco.Eta(), sv_reco.Phi(), dv_3d.Eta(), dv_3d.Phi()));
                      tmp_2mutrk_cand.sv_d3Dsig.push_back(distsv.distance(sv, pvv).significance());
                      tmp_2mutrk_cand.sv_ppdl3D.push_back(temp_sv_d3D*temp_sv_cosdphi_3D*sv.p4().M()/sv.p4().P());

                      Int_t temp_sv_nmu =0;

                      for(Vertex::trackRef_iterator itk = sv.tracks_begin(); itk != sv.tracks_end(); itk++) {
                        for(size_t imu = 0; imu < softmueta.size(); imu ++) {
                            if(deltaR(softmueta[imu], softmuphi[imu], (**itk).eta(), (**itk).phi())<0.01)
                              temp_sv_nmu ++;
                        }
                      }

                      tmp_2mutrk_cand.sv_nmu.push_back(temp_sv_nmu);
                      tmp_2mutrk_cand.sv_mass.push_back(sv.p4().M());
                      tmp_2mutrk_cand.sv_pt.push_back(sv.p4().Pt());
                      tmp_2mutrk_cand.sv_dz.push_back(abs(dz));
                      tmp_2mutrk_cand.sv_ntrk.push_back(sv.tracksSize());
                      tmp_2mutrk_cand.n_sv++;
                    }

                    tmp_2mutrk_cand.dzpv[0] = 1;
                    tmp_2mutrk_cand.dzpv[1] = 1;
                    tmp_2mutrk_cand.dzpv[2] = 1;

                    tmp_2mutrk_cand.dzpv[0] = abs(sortedMu[0].innerTrack()->dz(pv1P));
                    tmp_2mutrk_cand.dzpv[1] = abs(sortedMu[1].innerTrack()->dz(pv1P));
                    tmp_2mutrk_cand.dzpv[2] = abs(trk.dz(pv1P));

                    for(size_t jpv = 0; jpv < vertexes->size(); jpv++) {
                      if(jpv==ipv2)continue;
                      const Vertex & vi = (*vertexes)[jpv];
                      if(abs(sortedMu[0].innerTrack()->dz(vi.position()))<tmp_2mutrk_cand.mu_dz[0]) tmp_2mutrk_cand.dzpv[0]=-1;
                      if(abs(sortedMu[1].innerTrack()->dz(vi.position()))<tmp_2mutrk_cand.mu_dz[1]) tmp_2mutrk_cand.dzpv[1]=-1;
                      if(abs(sortedMu[2].innerTrack()->dz(vi.position()))<tmp_2mutrk_cand.mu_dz[2]) tmp_2mutrk_cand.dzpv[2]=-1;
                    }

                    ////////////////////
                    // Track Isolation
                    // How to decide if a track is associated with a certain PV ?
                    ////////////////////

                    double pttrk_tau = 0, pttrk_tau0p5 = 0,  pttrk_m1 = 0, pttrk_m2 = 0, pttrk_m3 = 0;

                    tmp_2mutrk_cand.mindca_iso = 999;
                    tmp_2mutrk_cand.mindca_iso05 = 999;

                    // Initialize track parameters
                    tmp_2mutrk_cand.ntrk_tau = 0; 
                    tmp_2mutrk_cand.ntrk_tau0p5 = 0; 
                    tmp_2mutrk_cand.ntrk_tau_b = 0; 
                    tmp_2mutrk_cand.ntrk_sum = 0;
                    for (size_t it = 0; it<3; ++it) tmp_2mutrk_cand.ntrk[it] = 0;  
                    tmp_2mutrk_cand.ntrk0p1 = 0; 
                    tmp_2mutrk_cand.ntrk0p2 = 0;
                    tmp_2mutrk_cand.ntrk0p5 = 0;
                    tmp_2mutrk_cand.maxdxy_pv0 = -1;

                    math::XYZPoint fvP = math::XYZPoint(fv.position().x(), fv.position().y(), fv.position().z());

                    for(size_t itk = 0; itk < tracks->size(); itk++) {
                      const reco::Track & t = (*tracks)[itk];
                      if(!(t.quality(TrackBase::tight)))continue;
                      if(deltaR(sortedMu[0].eta(), sortedMu[0].phi(), t.eta(), t.phi())<0.01)continue;
                      if(deltaR(sortedMu[1].eta(), sortedMu[1].phi(), t.eta(), t.phi())<0.01)continue;
                      if(deltaR(trk.eta(), trk.phi(), t.eta(), t.phi())<0.01)continue;

                      double dz = abs(t.dz(fvP));
                      double dxy = abs(t.dxy(fvP));
                      double dca_fv = sqrt(dz*dz+dxy*dxy);
                      double dr_tau = deltaR(t.eta(), t.phi(), vds.Eta(), vds.Phi());

                      // iso no. 1b - using pt_min, drtau_max of the 2 mu + trk
                      if(t.pt() > 0.33*tmp_2mutrk_cand.mu_pt_min && dr_tau < 3.*tmp_2mutrk_cand.mu_dsdR_max && dca_fv<0.05 ) {
                        pttrk_tau += t.pt();
                        tmp_2mutrk_cand.ntrk_tau++; // iso 3b
                        if(dca_fv<tmp_2mutrk_cand.mindca_iso) tmp_2mutrk_cand.mindca_iso=dca_fv; // iso 4b
                      } 

                      if(t.pt()<1.0) continue;  // was 1.2
                      // iso no. 1
                      if(dr_tau < 0.5 && dca_fv<0.05 ) {
                        pttrk_tau0p5 += t.pt();
                        tmp_2mutrk_cand.ntrk_tau0p5++; // iso 3
                        //if(dca_fv<mindca_iso05)mindca_iso05=dca_fv; // iso 4
                      }

                      if(dca_fv<0.05) tmp_2mutrk_cand.ntrk_tau_b++; // iso 3b
                      if(dca_fv<tmp_2mutrk_cand.mindca_iso05) tmp_2mutrk_cand.mindca_iso05=dca_fv; // iso 4

                      TransientTrack trkiso = theB->build(t);
                      ClosestApproachInRPhi cAppm1, cAppm2, cAppm3;
                      cAppm1.calculate(trkiso.initialFreeState(), t_trks[0].initialFreeState());
                      cAppm2.calculate(trkiso.initialFreeState(), t_trks[1].initialFreeState());
                      cAppm3.calculate(trkiso.initialFreeState(), t_trks[2].initialFreeState());
                      if(!(cAppm1.status()&&cAppm2.status()&&cAppm3.status())) continue;

                      // iso no. 2
                      if(deltaR(t.eta(), t.phi(), sortedMu[0].eta(), sortedMu[0].phi()) < 0.3 && cAppm1.distance() < 0.1) {// && dz1 < .3) 
                        tmp_2mutrk_cand.ntrk[0]++;
                        pttrk_m1 += t.pt();
                      }
                      if(deltaR(t.eta(), t.phi(), sortedMu[1].eta(), sortedMu[1].phi()) < 0.3 && cAppm2.distance() < 0.1) {//&& dz2 < .3) 
                        tmp_2mutrk_cand.ntrk[1]++;
                        pttrk_m2 += t.pt();
                      }
                      if(deltaR(t.eta(), t.phi(), trk.eta(), trk.phi()) < 0.3 && cAppm3.distance() < 0.1) {//&& dz3 < .3) 
                        tmp_2mutrk_cand.ntrk[2]++;
                        pttrk_m3 += t.pt();
                      }
                      if( (deltaR(t.eta(), t.phi(), sortedMu[0].eta(), sortedMu[0].phi()) < 0.3 && cAppm1.distance() < 0.1 )
                            ||(deltaR(t.eta(), t.phi(), sortedMu[1].eta(), sortedMu[1].phi()) < 0.3 && cAppm2.distance() < 0.1 )
                            ||(deltaR(t.eta(), t.phi(), trk.eta(), trk.phi()) < 0.3 && cAppm3.distance() < 0.1 )
                      ) tmp_2mutrk_cand.ntrk_sum++;


                      // displaced track counting
                      // only tracks consistent with PV
                      double dz_pv0=abs(t.dz(pv1P));
                      if(!(dz_pv0 < 1))continue;
                      double dxy_pv0 = abs(t.dxy(pv1P));
                      if(dxy_pv0>0.1) tmp_2mutrk_cand.ntrk0p1++;
                      if(dxy_pv0>0.2) tmp_2mutrk_cand.ntrk0p2++;
                      if(dxy_pv0>0.5) tmp_2mutrk_cand.ntrk0p5++;
                      if(dxy_pv0>tmp_2mutrk_cand.maxdxy_pv0) tmp_2mutrk_cand.maxdxy_pv0 = dxy_pv0;

                      tmp_2mutrk_cand.trkrel_tau = pttrk_tau/vds.Pt();
                      tmp_2mutrk_cand.trkrel_tau0p5 = pttrk_tau0p5/vds.Pt();
                      tmp_2mutrk_cand.trkrel[0] = pttrk_m1/sortedMu[0].pt(); 
                      tmp_2mutrk_cand.trkrel[1] = pttrk_m2/sortedMu[1].pt();
                      tmp_2mutrk_cand.trkrel[2] = pttrk_m3/trk.pt();
                      tmp_2mutrk_cand.trkrel_max = TMath::Max(tmp_2mutrk_cand.trkrel[0], TMath::Max(tmp_2mutrk_cand.trkrel[1], tmp_2mutrk_cand.trkrel[2]));
                    }
                }
                event_.diMuonTrk_cand_coll.push_back(tmp_2mutrk_cand);
                n2muTrk++;
              }

            }
        }

      }
    }
    event_.n3mu = n3mu;
    event_.n2muTrk = n2muTrk;
}

TransientVertex fitTriMuVertex(const edm::EventSetup & iSetup, TrackRef _trk1, TrackRef _trk2, TrackRef _trk3){
    ////////////////////
    // Fit 3mu vertex
    ////////////////////
    vector<TransientTrack> t_trks;
    KalmanVertexFitter kvf(true);
    TransientVertex fv = kvf.vertex(t_trks);
    ESHandle<TransientTrackBuilder> theB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
    t_trks.push_back(theB->build(_trk1));
    t_trks.push_back(theB->build(_trk2));
    t_trks.push_back(theB->build(_trk3));
    if(!fv.isValid()) { cout<<"[Tau23MuNtupleMaker]: Vertex Fit invalid!"<<endl; }
    return fv;
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Tau23MuNtupleMaker);
