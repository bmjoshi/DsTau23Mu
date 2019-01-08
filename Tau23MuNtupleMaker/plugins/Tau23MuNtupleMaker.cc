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

#include "Tau23MuNtupleMaker/Tools/src/Tau23MuTree.h"
#include "Utils.h"
#include "3mu_cand.h"
#include "pdg_info.h"

//#include "RecoMuon/MuonIdentification/src/MuonKinkFinder.cc"
//#include "TrackingTools/TrackRefitter/interface/RefitDirection.h"
//#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
//#include "TrackingTools/TrackFitters/interface/TrajectoryFitterRecord.h"
//#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
//#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
//#include "scnew.h"

using namespace reco;
using namespace l1t;
using namespace edm;
using namespace std;

#define dz_mumu_cut 0.5
#define dR_mumu_cut 0.8

class Tau23MuNtupleMaker : public edm::EDAnalyzer 
{
    public:

      Tau23MuNtupleMaker(const edm::ParameterSet &);

      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void beginRun(const edm::Run&, const edm::EventSetup&);
      virtual void beginJob();
      virtual void endJob();

    private:

      void fillAnalysis(const edm::Handle<reco::GenParticleCollection> & genParticles,
            const edm::Handle<std::vector<PileupSummaryInfo> > & puInfo,
            const edm::Handle<GenEventInfoProduct> & gen,
            const edm::Handle<edm::View<reco::Muon> > & muons,
            const edm::Handle<std::vector<reco::Track> >& tracks,
            const edm::Handle<std::vector<reco::Vertex> > & vertexes,
            const edm::Handle<reco::BeamSpot> & beamSpot,
            const edm::Handle<edm::TriggerResults> & triggerResults, 
            const edm::Handle<trigger::TriggerEvent> & triggerEvent,
            const edm::Handle<l1t::MuonBxCollection> & l1MuonBxColl,
            const edm::TriggerNames & triggerNames);
      void fillGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > &,
            const  edm::Handle<GenEventInfoProduct> &);

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


      // returns false in case the match is for a RPC chamber
      bool getMuonChamberId(DetId & id, muon_pog::MuonDetType & det, Int_t & r, Int_t & phi, Int_t & eta) const ;

      edm::EDGetTokenT<edm::TriggerResults> trigResultsToken_;
      edm::EDGetTokenT<trigger::TriggerEvent> trigSummaryToken_;

      std::string trigFilterCut_;
      std::string trigPathCut_;

      edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
      edm::EDGetTokenT<std::vector<reco::Track> > trackToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > primaryVertexToken_;
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

};


Tau23MuNtupleMaker::Tau23MuNtupleMaker( const edm::ParameterSet & cfg )
{
    // Branch flags

    edm::InputTag _getMu = cfg.getUntrackedParameter<>;


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

    m_minMuPtCut = cfg.getUntrackedParameter<double>("MinMuPtCut", 0.);
    m_minNMuCut  = cfg.getUntrackedParameter<int>("MinNMuCut",  0.);

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


void Tau23MuNtupleMaker::analyze (const edm::Event & ev, const edm::EventSetup &)
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

    event_.mets.pfMet   = -999; 
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


    // Fill general information
    // run, luminosity block, event
    event_.runNumber = ev.id().run();
    event_.luminosityBlockNumber = ev.id().luminosityBlock();
    event_.eventNumber = ev.id().event();

    eventId_.runNumber = ev.id().run();
    eventId_.luminosityBlockNumber = ev.id().luminosityBlock();
    eventId_.eventNumber = ev.id().event();

    // Fill GEN pile up information
    if (!ev.isRealData()) 
    {
      if (!pileUpInfoToken_.isUninitialized() &&
            !genInfoToken_.isUninitialized()) 
      {
        edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
        edm::Handle<GenEventInfoProduct> genInfo;

        if (ev.getByToken(pileUpInfoToken_, puInfo) &&
              ev.getByToken(genInfoToken_, genInfo) ) 
            fillGenInfo(puInfo,genInfo);
        else 
            edm::LogError("") << "[Tau23MuNtupleMaker]: Pile-Up Info collection does not exist !!!";
      }      
    }


    // Fill GEN particles information
    if (!ev.isRealData()) 
    {
      if (!genToken_.isUninitialized() ) 
      { 
        edm::Handle<reco::GenParticleCollection> genParticles;
        if (ev.getByToken(genToken_, genParticles)) 
            fillGenParticles(genParticles);
        else 
            edm::LogError("") << ">>> GEN collection does not exist !!!";
      }
    }

    if (ev.isRealData()) 
    {

      event_.bxId  = ev.bunchCrossing();
      event_.orbit = ev.orbitNumber();

      if (!scalersToken_.isUninitialized()) 
      { 
        edm::Handle<LumiScalersCollection> lumiScalers;
        if (ev.getByToken(scalersToken_, lumiScalers) && 
              lumiScalers->size() > 0 ) 
            event_.instLumi  = lumiScalers->begin()->instantLumi();
        else 
            edm::LogError("") << ">>> Scaler collection does not exist !!!";
      }
    }

    // Fill trigger information
    if (!trigResultsToken_.isUninitialized() &&
        !trigSummaryToken_.isUninitialized()) 
    {

      edm::Handle<edm::TriggerResults> triggerResults;
      edm::Handle<trigger::TriggerEvent> triggerEvent;

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

    // Get beam spot for muons
    edm::Handle<reco::BeamSpot> beamSpot;
    if (!beamSpotToken_.isUninitialized() ) 
    { 
      if (!ev.getByToken(beamSpotToken_, beamSpot)) 
        edm::LogError("") << "[Tau23MuNtupleMaker]: Beam spot collection not found !!!";
    }

    // Fill (raw) MET information: PF, PF charged, Calo    
    edm::Handle<reco::PFMETCollection> pfMet; 
    if(!pfMetToken_.isUninitialized()) 
    { 
      if (!ev.getByToken(pfMetToken_, pfMet)) 
        edm::LogError("") << "[Tau23MuNtupleMaker] PFMet collection does not exist !!!"; 
      else { 
        const reco::PFMET &iPfMet = (*pfMet)[0]; 
        event_.mets.pfMet = iPfMet.et(); 
      } 
    } 

    edm::Handle<reco::PFMETCollection> pfChMet; 
    if(!pfChMetToken_.isUninitialized()) 
    { 
      if (!ev.getByToken(pfChMetToken_, pfChMet)) 
        edm::LogError("") << "[Tau23MuNtupleMaker] PFChMet collection does not exist !!!"; 
      else { 
        const reco::PFMET &iPfChMet = (*pfChMet)[0]; 
        event_.mets.pfChMet = iPfChMet.et(); 
      } 
    } 

    edm::Handle<reco::CaloMETCollection> caloMet; 
    if(!caloMetToken_.isUninitialized()) 
    { 
      if (!ev.getByToken(caloMetToken_, caloMet)) 
        edm::LogError("") << "[Tau23MuNtupleMaker] CaloMet collection does not exist !!!"; 
      else { 
        const reco::CaloMET &iCaloMet = (*caloMet)[0]; 
        event_.mets.caloMet = iCaloMet.et(); 
      } 
    } 

    // Get muons  
    edm::Handle<edm::View<reco::Muon> > muons;
    if (!muonToken_.isUninitialized() ) 
    { 
      if (!ev.getByToken(muonToken_, muons)) 
        edm::LogError("") << "[Tau23MuNtupleMaker] Muon collection does not exist !!!";
    }

    //Get tracks
    edm::Handle<std::vector<reco::Track> > tracks;
    if (!trackToken_.isUninitialized())
    {
      if (!ev.getByToken(trackToken_, tracks))
        edm::LogError("") << "[Tau23MuNtupleMaker] Track collection does not exist !!!";
      else fillTracks(tracks,vertexes,beamSpot);

    }

    Int_t nGoodMuons = 0;
    eventId_.maxPTs.clear();
    // Fill muon information
    if (muons.isValid() && vertexes.isValid() && beamSpot.isValid()) 
    {
      nGoodMuons = fillMuons(muons,vertexes,beamSpot);
    }
    eventId_.nMuons = nGoodMuons;


    //Fill L1 informations
    edm::Handle<l1t::MuonBxCollection> l1s;
    if (!l1Token_.isUninitialized() )
    {
      if (!ev.getByToken(l1Token_, l1s))
        edm::LogError("") << "[Tau23MuNtupleMaker] L1 muon bx collection does not exist !!!";
      else {
        fillL1(l1s);
      }
    }

    if (nGoodMuons >= m_minNMuCut)
      tree_["tau23muTree"]->Fill();

}

void Tau23MuNtupleMaker::fillGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > & puInfo,
      const edm::Handle<GenEventInfoProduct> & gen)
{

    muon_pog::GenInfo genInfo;

    genInfo.trueNumberOfInteractions     = -1.;
    genInfo.actualNumberOfInteractions   = -1.;
    genInfo.genWeight = gen->weight() ;

    std::vector<PileupSummaryInfo>::const_iterator puInfoIt  = puInfo->begin();
    std::vector<PileupSummaryInfo>::const_iterator puInfoEnd = puInfo->end();

    for(; puInfoIt != puInfoEnd; ++puInfoIt) 
    {
      int bx = puInfoIt->getBunchCrossing();

      if(bx == 0) 
      { 
        genInfo.trueNumberOfInteractions   = puInfoIt->getTrueNumInteractions();
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

      muon_pog::GenParticle gensel;
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

            muon_pog::HLTObject hltObj;

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

        muon_pog::L1Muon l1part;
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
      const edm::Handle<reco::BeamSpot> & beamSpot)
{

    edm::View<reco::Muon>::const_iterator muonIt  = muons->begin();
    edm::View<reco::Muon>::const_iterator muonEnd = muons->end();

    for (; muonIt != muonEnd; ++muonIt) 
    {


      const reco::Muon& mu = (*muonIt);

      bool isGlobal      = mu.isGlobalMuon();
      bool isTracker     = mu.isTrackerMuon();
      bool isTrackerArb  = muon::isGoodMuon(mu, muon::TrackerMuonArbitrated); 
      bool isRPC         = mu.isRPCMuon();
      bool isStandAlone  = mu.isStandAloneMuon();
      bool isPF          = mu.isPFMuon();

      bool hasInnerTrack = !mu.innerTrack().isNull();
      bool hasTunePTrack = !mu.tunePMuonBestTrack().isNull();
      bool hasPickyTrack = !mu.pickyTrack().isNull();
      bool hasDytTrack = !mu.dytTrack().isNull();
      bool hasTpfmsTrack = !mu.tpfmsTrack().isNull();

      muon_pog::Muon ntupleMu;

      ntupleMu.pt     = mu.pt();
      ntupleMu.eta    = mu.eta();
      ntupleMu.phi    = mu.phi();
      ntupleMu.charge = mu.charge();
      ntupleMu.vx = mu.vx();
      ntupleMu.vy = mu.vy();
      ntupleMu.vz = mu.vz();

      ntupleMu.fits.push_back(muon_pog::MuonFit(mu.pt(),mu.eta(),mu.phi(),
              mu.charge(),mu.muonBestTrack()->ptError()));

      ntupleMu.fits.push_back(muon_pog::MuonFit(hasInnerTrack ? mu.innerTrack()->pt()  : -1000.,
              hasInnerTrack ? mu.innerTrack()->eta() : -1000.,
              hasInnerTrack ? mu.innerTrack()->phi() : -1000.,
              hasInnerTrack ? mu.innerTrack()->charge()  : -1000.,
              hasInnerTrack ? mu.innerTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(isStandAlone ? mu.outerTrack()->pt()  : -1000.,
              isStandAlone ? mu.outerTrack()->eta() : -1000.,
              isStandAlone ? mu.outerTrack()->phi() : -1000.,
              isStandAlone ? mu.outerTrack()->charge()  : -1000.,
              isStandAlone ? mu.outerTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(isGlobal ? mu.globalTrack()->pt()  : -1000.,
              isGlobal ? mu.globalTrack()->eta() : -1000.,
              isGlobal ? mu.globalTrack()->phi() : -1000.,
              isGlobal ? mu.globalTrack()->charge()  : -1000.,
              isGlobal ? mu.globalTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(hasTunePTrack ? mu.tunePMuonBestTrack()->pt()  : -1000.,
              hasTunePTrack ? mu.tunePMuonBestTrack()->eta() : -1000.,
              hasTunePTrack ? mu.tunePMuonBestTrack()->phi() : -1000.,
              hasTunePTrack ? mu.tunePMuonBestTrack()->charge()  : -1000.,
              hasTunePTrack ? mu.tunePMuonBestTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(hasPickyTrack ? mu.pickyTrack()->pt()  : -1000.,
              hasPickyTrack ? mu.pickyTrack()->eta() : -1000.,
              hasPickyTrack ? mu.pickyTrack()->phi() : -1000.,
              hasPickyTrack ? mu.pickyTrack()->charge()  : -1000.,
              hasPickyTrack ? mu.pickyTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(hasDytTrack ? mu.dytTrack()->pt()  : -1000.,
              hasDytTrack ? mu.dytTrack()->eta() : -1000.,
              hasDytTrack ? mu.dytTrack()->phi() : -1000.,
              hasDytTrack ? mu.dytTrack()->charge()  : -1000.,
              hasDytTrack ? mu.dytTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(hasTpfmsTrack ? mu.tpfmsTrack()->pt()  : -1000.,
              hasTpfmsTrack ? mu.tpfmsTrack()->eta() : -1000.,
              hasTpfmsTrack ? mu.tpfmsTrack()->phi() : -1000.,
              hasTpfmsTrack ? mu.tpfmsTrack()->charge()  : -1000.,
              hasTpfmsTrack ? mu.tpfmsTrack()->ptError() : -1000.));

      // Detector Based Isolation
      reco::MuonIsolation detIso03 = mu.isolationR03();

      ntupleMu.trackerIso = detIso03.sumPt;
      ntupleMu.EMCalIso   = detIso03.emEt;
      ntupleMu.HCalIso    = detIso03.hadEt;

      // PF Isolation
      reco::MuonPFIsolation pfIso04 = mu.pfIsolationR04();
      reco::MuonPFIsolation pfIso03 = mu.pfIsolationR03();

      ntupleMu.chargedHadronIso   = pfIso04.sumChargedHadronPt;
      ntupleMu.chargedHadronIsoPU = pfIso04.sumPUPt; 
      ntupleMu.neutralHadronIso   = pfIso04.sumNeutralHadronEt;
      ntupleMu.photonIso          = pfIso04.sumPhotonEt;

      ntupleMu.isGlobal     = isGlobal ? 1 : 0; 
      ntupleMu.isTracker    = isTracker ? 1 : 0; 
      ntupleMu.isTrackerArb = isTrackerArb ? 1 : 0; 
      ntupleMu.isRPC        = isRPC ? 1 : 0;
      ntupleMu.isStandAlone = isStandAlone ? 1 : 0;
      ntupleMu.isPF         = isPF ? 1 : 0;

      ntupleMu.nHitsGlobal     = isGlobal     ? mu.globalTrack()->numberOfValidHits() : -999; 
      ntupleMu.nHitsTracker    = isTracker    ? mu.innerTrack()->numberOfValidHits()  : -999; 
      ntupleMu.nHitsStandAlone = isStandAlone ? mu.outerTrack()->numberOfValidHits()  : -999;

      ntupleMu.glbNormChi2              = isGlobal      ? mu.globalTrack()->normalizedChi2() : -999; 
      ntupleMu.trkNormChi2         = hasInnerTrack ? mu.innerTrack()->normalizedChi2()  : -999; 
      ntupleMu.trkMuonMatchedStations   = isTracker     ? mu.numberOfMatchedStations()       : -999; 
      ntupleMu.glbMuonValidHits         = isGlobal      ? mu.globalTrack()->hitPattern().numberOfValidMuonHits()       : -999; 
      ntupleMu.trkPixelValidHits = hasInnerTrack ? mu.innerTrack()->hitPattern().numberOfValidPixelHits()       : -999; 
      ntupleMu.trkPixelLayersWithMeas   = hasInnerTrack ? mu.innerTrack()->hitPattern().pixelLayersWithMeasurement()   : -999; 
      ntupleMu.trkTrackerLayersWithMeas = hasInnerTrack ? mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() : -999; 

      ntupleMu.bestMuPtErr              = mu.muonBestTrack()->ptError(); 

      ntupleMu.trkValidHitFrac = hasInnerTrack           ? mu.innerTrack()->validFraction()       : -999; 
      ntupleMu.trkStaChi2      = isGlobal                ? mu.combinedQuality().chi2LocalPosition : -999; 
      ntupleMu.trkKink         = isGlobal                ? mu.combinedQuality().trkKink           : -999; 
      ntupleMu.muSegmComp      = (isGlobal || isTracker) ? muon::segmentCompatibility(mu)         : -999; 

      ntupleMu.isTrkMuOST               = muon::isGoodMuon(mu, muon::TMOneStationTight) ? 1 : 0; 
      ntupleMu.isTrkHP                  = hasInnerTrack && mu.innerTrack()->quality(reco::TrackBase::highPurity) ? 1 : 0; 

      if ( mu.isMatchesValid() && ntupleMu.isTrackerArb )
      {
        for ( reco::MuonChamberMatch match : mu.matches() )
        {
            muon_pog::ChambMatch ntupleMatch;

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
      ntupleMu.dzBest   = -999; 
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

      ntupleMu.isSoft    = 0;   
      ntupleMu.isTight   = 0;   
      ntupleMu.isHighPt  = 0;
      ntupleMu.isLoose   = muon::isLooseMuon(mu)  ? 1 : 0;   
      ntupleMu.isMedium  = muon::isMediumMuon(mu) ? 1 : 0;   

      if (vertexes->size() > 0)
      {
        const reco::Vertex & vertex = vertexes->at(0);

        dxy = isGlobal ? mu.globalTrack()->dxy(vertex.position()) :
            hasInnerTrack ? mu.innerTrack()->dxy(vertex.position()) : -1000;
        dz = isGlobal ? mu.globalTrack()->dz(vertex.position()) :
            hasInnerTrack ? mu.innerTrack()->dz(vertex.position()) : -1000;

        ntupleMu.dxyBest  = mu.muonBestTrack()->dxy(vertex.position()); 
        ntupleMu.dzBest   = mu.muonBestTrack()->dz(vertex.position()); 
        if(hasInnerTrack) { 
            ntupleMu.dxyInner = mu.innerTrack()->dxy(vertex.position()); 
            ntupleMu.dzInner  = mu.innerTrack()->dz(vertex.position()); 
        } 

        ntupleMu.isSoft    = muon::isSoftMuon(mu,vertex)   ? 1 : 0;   
        ntupleMu.isTight   = muon::isTightMuon(mu,vertex)  ? 1 : 0;   
        ntupleMu.isHighPt  = muon::isHighPtMuon(mu,vertex) ? 1 : 0;

      }

      ntupleMu.dxy    = dxy;
      ntupleMu.dz     = dz;
      ntupleMu.edxy   = isGlobal ? mu.globalTrack()->dxyError() : hasInnerTrack ? mu.innerTrack()->dxyError() : -1000;
      ntupleMu.edz    = isGlobal ? mu.globalTrack()->dzError()  : hasInnerTrack ? mu.innerTrack()->dzError() : -1000;

      ntupleMu.dxybs  = dxybs;
      ntupleMu.dzbs   = dzbs;

      if(mu.isTimeValid()) { 
        ntupleMu.muonTimeDof = mu.time().nDof; 
        ntupleMu.muonTime    = mu.time().timeAtIpInOut; 
        ntupleMu.muonTimeErr = mu.time().timeAtIpInOutErr; 
      } 
      else { 
        ntupleMu.muonTimeDof = -999; 
        ntupleMu.muonTime    = -999; 
        ntupleMu.muonTimeErr = -999; 
      } 

      if(mu.rpcTime().nDof > 0) { 
        ntupleMu.muonRpcTimeDof = mu.rpcTime().nDof; 
        ntupleMu.muonRpcTime    = mu.rpcTime().timeAtIpInOut; 
        ntupleMu.muonRpcTimeErr = mu.rpcTime().timeAtIpInOutErr; 
      } 
      else { 
        ntupleMu.muonRpcTimeDof = -999; 
        ntupleMu.muonRpcTime    = -999; 
        ntupleMu.muonRpcTimeErr = -999; 
      } 

      // asking for a TRK or GLB muon with minimal pT cut
      // ignoring STA muons in this logic
      if ( m_minMuPtCut < 0 ||
            (
           (isTracker || isGlobal || isStandAlone) &&
           (ntupleMu.fitPt(muon_pog::MuonFitType::DEFAULT) > m_minMuPtCut ||
            ntupleMu.fitPt(muon_pog::MuonFitType::GLB)     > m_minMuPtCut ||
            ntupleMu.fitPt(muon_pog::MuonFitType::TUNEP)   > m_minMuPtCut ||
            ntupleMu.fitPt(muon_pog::MuonFitType::INNER)   > m_minMuPtCut ||
            ntupleMu.fitPt(muon_pog::MuonFitType::STA)     > m_minMuPtCut ||
            ntupleMu.fitPt(muon_pog::MuonFitType::PICKY)   > m_minMuPtCut ||
            ntupleMu.fitPt(muon_pog::MuonFitType::DYT)     > m_minMuPtCut ||
            ntupleMu.fitPt(muon_pog::MuonFitType::TPFMS)   > m_minMuPtCut)
            )
       )
      {
        event_.muons.push_back(ntupleMu);

        std::vector<Float_t> PTs = {ntupleMu.fitPt(muon_pog::MuonFitType::DEFAULT),
            ntupleMu.fitPt(muon_pog::MuonFitType::GLB),
            ntupleMu.fitPt(muon_pog::MuonFitType::TUNEP),
            ntupleMu.fitPt(muon_pog::MuonFitType::INNER),
            ntupleMu.fitPt(muon_pog::MuonFitType::PICKY),
            ntupleMu.fitPt(muon_pog::MuonFitType::DYT),
            ntupleMu.fitPt(muon_pog::MuonFitType::TPFMS)};
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

      muon_pog::Track ntupleTr;

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

bool Tau23MuNtupleMaker::getMuonChamberId(DetId & id, muon_pog::MuonDetType & det,
      Int_t & r, Int_t & phi, Int_t & eta) const
{

    if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT)
    {
      DTChamberId dtId(id.rawId());  

      det = muon_pog::MuonDetType::DT;
      r   = dtId.station();
      phi = dtId.sector();
      eta = dtId.wheel();

      return true;
    }

    if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC)
    {
      CSCDetId cscId(id.rawId());

      det = muon_pog::MuonDetType::CSC;
      r   = cscId.station() * cscId.zendcap();
      phi = cscId.chamber();
      eta = cscId.ring();

      return true;
    }

    return false;

}

void Tau23MuNtupleMaker::fillAnalysisTree(const edm::Handle<reco::GenParticleCollection> & genParticles,
      const edm::Handle<std::vector<PileupSummaryInfo> > & puInfo,
      const edm::Handle<GenEventInfoProduct> & gen,
      const edm::Handle<edm::View<reco::Muon> > & muons,
      const edm::Handle<std::vector<reco::Track> >& tracks,
      const edm::Handle<std::vector<reco::Vertex> > & vertexes,
      const edm::Handle<reco::BeamSpot> & beamSpot,
      const edm::Handle<edm::TriggerResults> & triggerResults, 
      const edm::Handle<trigger::TriggerEvent> & triggerEvent,
      const edm::Handle<l1t::MuonBxCollection> & l1MuonBxColl,
      const edm::TriggerNames & triggerNames){

    tau23mu::3mu_cand tmp_cand;
    tau23mu::2muTrk_cand tmp_2mutrk_cand;

    reco::Track::const_iterator trackIt = tracks->begin();
    reco::Track::const_iterator trackEnd = tracks->end();

    edm::View<reco::Muon>::const_iterator muonIt  = muons->begin();
    edm::View<reco::Muon>::const_iterator muonEnd = muons->end();

    Int_t mu_idx = -1;
    Int_t tr_idx = -1;

    for (; muonIt != muonEnd; ++muonIt) 
    {
      const reco::Muon& mu = (*muonIt);
      mu_idx++;
      if (mu.pt()<2 || abs(mu.eta())>2.4) continue;
      if (mu.isPFMuon() || mu.isGlobalMuon()) goodMuIndex.push_back(mu_idx); // Store indices of all good muons
    }

    for (; trackIt != trackEnd; ++trackIt){
      const reco::Track& trk = (*trackIt);
      trk_idx++;

      for (; muonIt != muonEnd; ++muonIt){
        const reco::Muon& mu = (*muonIt);
        if (abs(trk.eta()-mu.eta())<0.001 && abs(trk.phi()-mu.phi())<0.001) muMatchedTrack.push_back(trk_idx);
      }
    }


    if (goodMuIndex.size()<(do2mu_?2:3)) return; // Skip, if the number of good muons is less than 3(2) for Tau3Mu(DsPhiPi) analysis

    for (size_t i=0; i<goodMuIndex.size()-1; ++i){
      const reco::Muon & mu1 = (*muons)[goodMuIndex[i]]; // Select first good muon

      for (size_t j=i+1; j<goodMuIndex.size(); ++j){
        const reco::Muon & mu2 = (*muons)[goodMuonIndex[j]]; // Select second good muon

        double dz_mu1mu2 = abs(mu1.InnerTracker()->dz(beamSpotHandle->position())-mu2.InnerTracker()->dz(beamSpotHandle->Position()));
        if ( dz_mu1mu2 > dz_mu1mu2_cut ) continue; // dz cut 

        double dR_mu1mu2 = deltaR(mu1.eta(), mu2.eta(), mu1.phi(), mu2.phi());
        if ( dR_mu1mu2 > dR_mumu_cut ) continue; // dR cut

        if ( goodMuonIndex.size() < j-1){

            for (size_t k=j+1; k<goodMuonIndex.size(); ++k){

              condt reco::Muon & mu3 = (*muons)[goodMuonIndex[k]];

              size_t nPt2p5 = 0;
              if (mu1.pt()>2.5) nPt2p5++;
              if (mu2.pt()>2.5) nPt2p5++;
              if (mu3.pt()>2.5) nPt2p5++;

              if (nPt2p5<2) continue; // Require at least 2 muons with Pt more than 2.5 GeV

              double dz_mu2mu3 = abs(mu2.InnerTracker()->dz(beamSpotHandle->position())-mu3.InnerTracker()->dz(beamSpotHandle->Position()));
              double dz_mu1mu3 = abs(mu1.InnerTracker()->dz(beamSpotHandle->position())-mu3.InnerTracker()->dz(beamSpotHandle->Position()));

              if (dz_mu2mu3 > dz_mumu_cut || dz_mu1mu3 > dz_mumu_cut) continue; // dz cut

              double dR_mu1mu2 = deltaR(mu1.eta(), mu2.eta(), mu1.phi(), mu2.phi()); 
              double dR_mu1mu2 = deltaR(mu1.eta(), mu2.eta(), mu1.phi(), mu2.phi());

              if( dR_mu2mu3 > dR_mumu_cut || dR_mu1mu3 > dR_mumu_cut ) continue; // dR cut

              if ( mu1.charge()+mu2.charge()+mu3.charge() > 1) continue; // Require sum of charges to be +/- 1


              // Build tracks  and re-fit the vertex
              TrackRef trk1 = mu1.innerTrack();
              TrackRef trk2 = mu2.innerTrack();
              TrackRef trk2 = mu3.innerTrack();

              auto fv = tau23mu::fitTriMuVertex(trk1, trk2, trk3);
              double fvnC2_tmp = fv.totalChiSquared()/fv.degreesOfFreedom();

              if (fvnC2_tmp > 3muFitVtx_nC2_cut ) continue; // Store all the good candidates

              // Sort the muons by P
              vector<reco::Muons> sortedMu;

              sortedMu.push_back(maxPMu(maxPMu(mu1,mu2),mu3));
              sortedMu.push_back(minPMu(minPMu(maxPMu(mu1,mu2),maxPMu(mu2,mu3)),maxPMu(mu3,mu1)));
              sortedMu.push_back(minPMu(minPMu(mu1,mu2),mu3));

              // Fill tree with trimuon information
              for (size_t it=0; it<3; ++it){
                tmp_cand.mu_p[it] = sortedMu[it].p();
                tmp_cand.mu_pt[it] = sotredMu[it].pt();
                tmp_cand.mu_pterr[it] = sortedMu[it].muonBestTrack()->ptError();
                tmp_cand.mu_bestPt[it] = sortedMu[it].muonBestTrack->pt();
                tmp_cand.mu_eta[it] = sotredMu[it].eta();
                tmp_cand.mu_phi[it] = sotredMu[it].phi();
                tmp_cand.mu_charge[it] = sortedMu[it].charge();
                tmp_cand.mu_vx[it] = sortedMu[it].vx();
                tmp_cand.mu_vy[it] = sortedMu[it].vy();
                tmp_cand.mu_vz[it] = sortedMu[it].vz();

                // muon ID
                bool _isGlobal = sortedMu[it].isGlobalMuon();
                bool _isTracker = sortedMu[it].isTrackerMuon();
                bool _isTrackerArb = muon::isGoodMuon(sortedMu[it].TrackerArbitrated);
                bool _isRPC = sortedMu[it].isRPCMuon();
                bool _isStandAlone = sortedMu[it].isStandAloneMuon();
                bool _isPF = sortedMu[it].isPFMuon();

                bool _hasInnerTrack = !mu.innerTrack().isNull();
                bool _hasTunePTrack = !mu.tunePMuonBestTrack().isNull();
                bool _hasPickyTrack = !mu.pickyTrack().isNull();
                bool _hasDytTrack = !mu.dytTrack().isNull();
                bool _hasTpfmsTrack = !mu.tpfmsTrack().isNull();

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

                tmp_cand.mu_uSta[it] = sortedMu[it].combinedQuality().updatedSta;
                tmp_cand.mu_trkKink[it] = sortedMu[it].combinedQuality().trkKink;
                tmp_cand.mu_log_glbKink[it] = TMath::Log(2+sortedMu[it].combinedQuality().glbKink);
                tmp_cand.mu_tRC2[it] = sortedMu[it].combinedQuality().trkRelChi2;
                tmp_cand.mu_sRC2[it] = sortedMu[it].combinedQuality().staRelChi2;
                tmp_cand.mu_Chi2LP[it] = sortedMu[it].combinedQuality().chi2LocalPosition;
                tmp_cand.mu_Chi2LM[it] = sortedMu[it].combinedQuality().chi2LocalMomentum;
                tmp_cand.mu_localDist[it] = sortedMu[it].combinedQuality().localDistance;
                tmp_cand.mu_glbDEP[it] = sortedMu[it].combinedQuality().globalDeltaEtaPhi;
                tmp_cand.mu_tightMatch[it] = sortedMu[it].combinedQuality().tightMatch;
                tmp_cand.mu_glbTrkProb[it] = sortedMu[it].combinedQuality.glbTrackProbability;
                tmp_cand.mu_calEM[it] = sortedMu[it].calEnergy().em;
                tmp_cand.mu_calEMS9[it] = sortedMu[it].calEnergy().emS9;
                tmp_cand.mu_calEMS25[it] = sortedMu[it].calEnergy().emS25;
                tmp_cand.mu_calHad[it] = sortedMu[it].calEnergy().had;
                tmp_cand.mu_calHadS9[it] = sortedMu[it].calEnergy.hadS9;;
                tmp_cand.mu_nOMS[it] = sortedMu[it].numberOfMatchedStations();
                tmp_cand.mu_nOM[it] = sortedMu[it].numberOfMatches(reco::Muon::segmentArbitraction);
                tmp_cand.mu_comp2d[it] = muon::isGoodMuon(sortedMu[it], muon::TM2DCompatibilityTight);
                tmp_cand.mu_calocomp[it] = muon::caloCompatibility(sortedMu[it]);
                tmp_cand.mu_segcomp[it] = muon::segmentCompatibility(sortedMu[it]);

                // Muon outer track info
                tmp_cand.mu_oTrk_p[it] = _isStandAlone ? sortedMu[it].outerTrack()->p():-999;
                tmp_cand.mu_oTrk_pt[it] = _isStandAlone ? sortedMu[it].outerTrack()->pt():-999;
                tmp_cand.mu_oTrk_eta[it] = _isStandAlone ? sortedMu[it].outerTrack()->eta():-999;
                tmp_cand.mu_oTrk_phi[it] = _isStandAlone ? sortedMu[it].outerTrack()->phi(): -999;
                tmp_cand.mu_oTrk_outernC2[it] = _isStandAlone ? sortedMu[it].outerTrack()->normalizedChi2():-999;
                tmp_cand.mu_oTrk_MSWVH_[it] = _isStandAlone ? sortedMu[it].outerTrack()->hitPattern().muonStationsWithValidHits():-999;
                tmp_cand.mu_oTrk_qprod[it] = _isStandAlone ? sortedMu[it].outerTrack()->charge()*sortedMu[it].innerTrack()->charge():-999;

                // Muon global track info
                tmp_cand.mu_glb_p[it] = _isGlobal ? sortedMu[it].globalTrack()->p():-999;
                tmp_cand.mu_glb_eta[it] = _isGlobal ? sortedMu[it].globalTrack()->eta():-999;
                tmp_cand.mu_glb_phi[it] = _isGlobal ? sortedMu[it].globalTrack()->phi():-999;
                tmp_cand.mu_glb_nC2[it] = _isGlobal ? sortedMu[it].globalTrack()->normalizedChi2():-999;
                tmp_cand.mu_glb_nOVMH[it] = _isGlobal ? sortedMu[it].globalTrack()->hitPattern().numberOfValidMuonHits():-999;

                // muon inner track information
                tmp_cand.mu_inTrk_p[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->p():-999;
                tmp_cand.mu_inTrk_pt[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->pt():-999;
                tmp_cand.mu_inTrk_eta[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->eta():-999;
                tmp_cand.mu_inTrk_phi[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->phi():-999;
                tmp_cand.mu_inTrk_nC2[i] = _hasInnerTrack ? sortedMu[i].innerTrack()->normalizedChi2():-999;
                tmp_cand.mu_inTrk_trkLayWithMeas[it] = _hasInnerTrack ? sortedMu[it].innerTracker()->hitPattern().trackerLayersWithMeasurement():-999;
                tmp_cand.mu_inTrk_pixLayWithMeas[it] = _hasInnerTrack ? sortedMu[it].innerTracker()->hitPattern().pixelLayersWithMeasurement():-999;
                tmp_cand.mu_inTrk_nValidTrkHits[it] = _hasInnerTrack ? sortedMu[it].innerTracker()->hitPattern().numberOfValidTrackerHits():-999;
                tmp_cand.mu_inTrk_nValidPixHits[it] = _hasInnerTrack ? sortedMu[it].innerTracker()->hitPattern().numberOfValidTrackerHits():-999;
                tmp_cand.mu_inTrk_validFraction[it] = _hasInnerTrack ? sortedMu[it].innerTracker()->validFraction():-999;
                tmp_cand.mu_inTrk_nLostTrkHits[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfLostTrackerHits(hitPattern::TRACK_HITS):-999;
                tmp_cand.mu_inTrk_nLostTrkHits_in[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfLostTrackerHits(hitPattern::MISSING_INNER_HITS):-999;
                tmp_cand.mu_inTrk_nLostTrkHits_out[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfLostTrackerHits(hitPattern::MISSING_OUTER_HITS):-999;
                tmp_cand.mu_inTrk_HP[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->quality(TrackBase::highPurity);

                // Isolation variables
                reco::MuonIsolation IsoR03 = sortedMu[it].isolationR03;

                tmp_cand.mu_IsoR03_sumPt[it] = IsoR03.sumPt;
                tmp_cand.mu_IsoR03_nTrks[it] = IsoR03.nTracks;
                tmp_cand.mu_IsoR03_emEt[it] = IsoR03.emEt;
                tmp_cand.mu_IsoR03_hadEt[it] = IsoR03.hadEt;
                tmp_cand.mu_IsoR03_emVetoEt[it] = IsoR03.emVetoEt;
                tmp_cand.mu_IsoR03_hadVetoEt[it] = IsoR03.hadVetoEt;

                // PF isolation variables
                reco::MuonPFIsolation pfIso04 = sortedMu[it].pfIsolationR04();
                reco::MuonPFIsolation pfIso03 = sortedMu[it].pfIsolationR03();

                tmp_cand.chargedHadronIso   = pfIso04.sumChargedHadronPt;
                tmp_cand.chargedHadronIsoPU = pfIso04.sumPUPt; 
                tmp_cand.neutralHadronIso   = pfIso04.sumNeutralHadronEt;
                tmp_cand.photonIso          = pfIso04.sumPhotonEt;

                tmp_cand.fits.push_back(tau23mu::MuonFit(sortedMu[it].pt(),sortedMu[it].eta(),sortedMu[it].phi(),
                        sortedMu[it].charge(),sortedMu[it].muonBestTrack()->ptError()));

                tmp_cand.fits.push_back(tau23mu::MuonFit(hasInnerTrack ? sortedMu[it].innerTrack()->pt()  : -1000.,
                        hasInnerTrack ? sortedMu[it].innerTrack()->eta() : -1000.,
                        hasInnerTrack ? sortedMu[it].innerTrack()->phi() : -1000.,
                        hasInnerTrack ? sortedMu[it].innerTrack()->charge()  : -1000.,
                        hasInnerTrack ? sortedMu[it].innerTrack()->ptError() : -1000.));

                tmp_cand.fits.push_back(tau23mu::MuonFit(isStandAlone ? sortedMu[it].outerTrack()->pt()  : -1000.,
                        isStandAlone ? sortedMu[it].outerTrack()->eta() : -1000.,
                        isStandAlone ? sortedMu[it].outerTrack()->phi() : -1000.,
                        isStandAlone ? sortedMu[it].outerTrack()->charge()  : -1000.,
                        isStandAlone ? sortedMu[it].outerTrack()->ptError() : -1000.));

                tmp_cand.fits.push_back(tau23mu::MuonFit(isGlobal ? sortedMu[it].globalTrack()->pt()  : -1000.,
                        isGlobal ? sortedMu[it].globalTrack()->eta() : -1000.,
                        isGlobal ? sortedMu[it].globalTrack()->phi() : -1000.,
                        isGlobal ? sortedMu[it].globalTrack()->charge()  : -1000.,
                        isGlobal ? sortedMu[it].globalTrack()->ptError() : -1000.));

                tmp_cand.fits.push_back(tau23mu::MuonFit(hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->pt()  : -1000.,
                        hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->eta() : -1000.,
                        hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->phi() : -1000.,
                        hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->charge()  : -1000.,
                        hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->ptError() : -1000.));

                tmp_cand.fits.push_back(tau23mu::MuonFit(hasPickyTrack ? sortedMu[it].pickyTrack()->pt()  : -1000.,
                        hasPickyTrack ? sortedMu[it].pickyTrack()->eta() : -1000.,
                        hasPickyTrack ? sortedMu[it].pickyTrack()->phi() : -1000.,
                        hasPickyTrack ? sortedMu[it].pickyTrack()->charge()  : -1000.,
                        hasPickyTrack ? sortedMu[it].pickyTrack()->ptError() : -1000.));

                tmp_cand.fits.push_back(tau23mu::MuonFit(hasDytTrack ? sortedMu[it].dytTrack()->pt()  : -1000.,
                        hasDytTrack ? sortedMu[it].dytTrack()->eta() : -1000.,
                        hasDytTrack ? sortedMu[it].dytTrack()->phi() : -1000.,
                        hasDytTrack ? sortedMu[it].dytTrack()->charge()  : -1000.,
                        hasDytTrack ? sortedMu[it].dytTrack()->ptError() : -1000.));

                tmp_cand.fits.push_back(tau23mu::MuonFit(hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->pt()  : -1000.,
                        hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->eta() : -1000.,
                        hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->phi() : -1000.,
                        hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->charge()  : -1000.,
                        hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->ptError() : -1000.));

                tmp_cand.mu_dxyBest  = -999; 
                tmp_cand.mu_dzBest   = -999; 
                tmp_cand.mu_dxyInner = -999; 
                tmp_cand.mu_dzInner  = -999; 

                tmp_cand.mu_dxybs = _isGlobal ? mu.globalTrack()->dxy(beamSpot->position()) :
                    _hasInnerTrack ? mu.innerTrack()->dxy(beamSpot->position()) : -1000;
                tmp_cand.mu_dzbs  = _isGlobal ? mu.globalTrack()->dz(beamSpot->position()) :
                    _hasInnerTrack ? mu.innerTrack()->dz(beamSpot->position()) : -1000;

                double dxy = -1000.;
                double dz  = -1000.;

                tmp_cand.isSoft    = 0;   
                tmp_cand.isTight   = 0;   
                tmp_cand.isHighPt  = 0;
                tmp_cand.isLoose   = muon::isLooseMuon(sortedMu[it])  ? 1 : 0;   
                tmp_cand.isMedium  = muon::isMediumMuon(sortedMu[it]) ? 1 : 0;   
                if ( mu.isMatchesValid() && _isTrackerArb )
                {
                    for ( reco::MuonChamberMatch match : mu.matches() )
                    {
                      muon_pog::ChambMatch ntupleMatch;

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

                        tmp_cand.matches.push_back(ntupleMatch);
                      }
                    }
                }
                tmp_cand.mu_isoPflow04 = (pfIso04.sumChargedHadronPt+ 
                      std::max(0.,pfIso04.sumPhotonEt+pfIso04.sumNeutralHadronEt - 0.5*pfIso04.sumPUPt)) / mu.pt();

                tmp_cand.mu_isoPflow03 = (pfIso03.sumChargedHadronPt+ 
                      std::max(0.,pfIso03.sumPhotonEt+pfIso03.sumNeutralHadronEt - 0.5*pfIso03.sumPUPt)) / mu.pt();

                if (vertexes->size() > 0) {
                    const reco::Vertex & vertex = vertexes->at(0);

                    tmp_cand.mu_dxy[it] = _isGlobal ? sortedMu[it].globalTrack()->dxy(vertex.position()) :
                      _hasInnerTrack ? sortedMu[it].innerTrack()->dxy(vertex.position()) : -1000;
                    tmp_cand.mu_dz[it] = _isGlobal ? sortedMu[it].globalTrack()->dz(vertex.position()) :
                      _hasInnerTrack ? sortedMu[it].innerTrack()->dz(vertex.position()) : -1000;

                    tmp_cand.mu_dxyBest[it]  = sortedMu[it].muonBestTrack()->dxy(vertex.position()); 
                    tmp_cand.mu_dzBest[it]   = sortedMu[it].muonBestTrack()->dz(vertex.position()); 
                    if(_hasInnerTrack) { 
                      tmp_cand.mu_dxyInner[it] = sortedMu[it].innerTrack()->dxy(vertex.position()); 
                      tmp_cand.mi_dzInner[it]  = sortedMu[it].innerTrack()->dz(vertex.position()); 
                    } 

                    tmp_cand.isSoft[it]    = muon::isSoftMuon(sortedMu[it],vertex)   ? 1 : 0;   
                    tmp_cand.isTight[it]   = muon::isTightMuon(sortedMu[it],vertex)  ? 1 : 0;   
                    tmp_cand.isHighPt[it]  = muon::isHighPtMuon(sortedMu[it],vertex) ? 1 : 0;

                }
                if(sortedMu[it].isTimeValid()) { 
                    tmp_cand.mu_TimeDof = sortedMu[it].time().nDof; 
                    tmp_cand.mu_Time    = sortedMu[it].time().timeAtIpInOut; 
                    tmp_cand.mu_TimeErr = sortedMu[it].time().timeAtIpInOutErr; 
                } 
                else { 
                    tmp_cand.muonTimeDof = -999; 
                    tmp_cand.muonTime    = -999; 
                    tmp_cand.muonTimeErr = -999; 
                } 

                if(sortedMu[it].rpcTime().nDof > 0) { 
                    tmp_cand.mu_RpcTimeDof[it] = sortedMu[it].rpcTime().nDof; 
                    tmp_cand.mu_RpcTime[it]    = sortedMu[it].rpcTime().timeAtIpInOut; 
                    tmp_cand.mu_RpcTimeErr[it] = sortedMu[it].rpcTime().timeAtIpInOutErr; 
                } 
                else { 
                    tmp_cand.mu_RpcTimeDof[it] = -999; 
                    tmp_cand.mu_RpcTime[it]   = -999; 
                    tmp_cand.mu_RpcTimeErr[it] = -999; 
                } 


                tmp_cand.mu_dxy[it]  = _isGlobal ? sortedMu[it].globalTrack()->dxy(vertex.position());
                tmp_cand.mu_dz[it]  = _isGlobal ? sortedMu[it].globalTrack()->dxy(vertex.position());
                tmp_cand.mu_edxy[it]   = _isGlobal ? sortedMu[it].globalTrack()->dxyError() : _hasInnerTrack ? sortedMu[it].innerTrack()->dxyError() : -1000;
                tmp_cand.mu_edz[it]    = _isGlobal ? sortedMu[it].globalTrack()->dzError()  : _hasInnerTrack ? sortedMu[it].innerTrack()->dzError() : -1000;

                tmp_cand.mu_TMOST[it] = muon::isGoodMuon(mu[i], muon::TMOneStationTight);
                tmp_cand.mu_TMOSAT[it] = muon::isGoodMuon(mu[i], muon::TMOneStationAngTight);
                tmp_cand.mu_TMLST[it] = muon::isGoodMuon(mu[i], muon::TMLastStationTight);
                tmp_cand.mu_TMLSAT[it] = muon::isGoodMuon(mu[i], muon::TMLastStationAngTight);
                tmp_cand.mu_TMLSOLPT[it] = muon::isGoodMuon(mu[i], muon::TMLastStationOptimizedLowPtTight);
                tmp_cand.mu_TMLSOBLPT[it] = muon::isGoodMuon(mu[i], muon::TMLastStationOptimizedBarrelLowPtTight);

                // Tau-Muon variables
                tmp_cand.mu_taudR[it] = deltaR(sortedMu[it].eta(), sortedMu[it].phi(), vtau.Eta(), vtau.Phi());
              }

              tmp_cand.mu_pt_max =  TMath::Max(TMath::Max(mu1.pt(),mu2.pt()),mu3.pt()); //return highest pt muon 
              tmp_cand.mu_eta_max = TMath::Max(TMath::Max(mu1.eta(),mu2.eta()),mu3.eta()); //return highest eta muon 
              tmp_cand.mu_pt_max =  TMath::Min(TMath::Min(mu1.pt(),mu2.pt()),mu3.pt()); //return lowest pt muon 
              tmp_cand.mu_eta_max = TMath::Min(TMath::Min(mu1.eta(),mu2.eta()),mu3.eta()); //return lowest eta muon
				  
				  tmp_cand.mu_taudR_max = TMath::Max(TMath::Max(tmp_cand.mu_taudR[0],tmp_cand.mu_taudR[1]),tmp_cand.mu_taudR[2]);
				  tmp_cand.mu_trkLayWithMeas_max = TMath::Max(TMath::Max(tmp_cand.mu_inTrk_trkLayWithMeas[0],tmp_cand.mu_inTrk_trkLayWithMeas[1]),tmp_cand.mu_inTrk_trkLayWithMeas[2]);
				  tmp_cand.mu_n

              // Re-assign and store dR and dz with sorted muons
              tmp_cand.dR_mu1mu2 = deltaR(sortedMu[0].eta(), sortedMu[1].eta(), sortedMu[0].phi(), sortedMu[1].phi());
              tmp_cand.dR_mu2mu3 = deltaR(sortedMu[1].eta(), sortedMu[2].eta(), sortedMu[1].phi(), sortedMu[2].phi());
              tmp_cand.dR_mu1mu3 = deltaR(sortedMu[0].eta(), sortedMu[2].eta(), sortedMu[0].phi(), sortedMu[2].phi());

              tmp_cand.dz_mu1mu2 = abs(sortedMu[0].InnerTracker()->dz(beamSpotHandle->position())-sortedMu[1].InnerTracker()->dz(beamSpotHandle->Position()));
              tmp_cand.dz_mu2mu3 = abs(sortedMu[1].InnerTracker()->dz(beamSpotHandle->position())-sortedMu[2].InnerTracker()->dz(beamSpotHandle->Position()));
              tmp_cand.dz_mu1mu3 = abs(sortedMu[0].InnerTracker()->dz(beamSpotHandle->position())-sortedMu[2].InnerTracker()->dz(beamSpotHandle->Position()));

              tmp_cand.mumu_dR_12 = dR_mu1mu2;
              tmp_cand.mumu_dR_23 = dR_mu2mu3;
              tmp_cand.mumu_dR_31 = dR_mu1mu3;

              tmp_cand.mumu_dz_12 = dz_mu1mu2;
              tmp_cand.mumu_dz_23 = dz_mu2mu3;
              tmp_cand.mumu_dz_31 = dz_mu1mu3;

              // Compute 3 mu invariant mass
              TLorentzVector vec_mu1, vec_mu2, vec_mu3;

              vec_mu1.SetPtEtaPhiM(sortedMu[0].pt(), sortedMu[0].eta(), sortedMu[0].phi(), MU_MASS);
              vec_mu2.SetPtEtaPhiM(sortedMu[1].pt(), sortedMu[1].eta(), sortedMu[1].phi(), MU_MASS);
              vec_mu3.SetPtEtaPhiM(sortedMu[2].pt(), sortedMu[2].eta(), sortedMu[2].phi(), MU_MASS);

              tmp_cand.M3mu.push_back((vec_mu1+vec_mu2+vec_mu3).M());
              tmp_cand.M_mu1mu2.push_back((vec_mu1+vec_mu2).M());
              tmp_cand.M_mu2mu3.push_back((vec_mu2+vec_mu3).M());
              tmp_cand.M_mu3mu1.push_back((vec_mu1+vec_mu3).M());

              tmp_cand.3muFitVtx_nC2.push_back(fvnC2_tmp);
              n_3mu++;
            }
        }
        //Muon variables from 3mu_cand

        std::vector<HitInfo> hits;
        // End of muon variables

        // Find 2mu+1track
        if (do2mu_ && mu1.pt()>2.5 && mu2.pt()>2.5){

            size_t trk_idx=-1, tmp_iter=0;

            // sort muons in the order of momentum
            reco::muon sortedMu[2];
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
                TrackRef trk3 = trk;

                TransientVertex fv = fitTriMuVertex(trk1, trk2, trk3);
                if(!fv.isValid()) continue;
                double fv_tC = fv.totalChiSquared();
                double fv_dOF = fv.degreesOfFreedom();
                double fv_nC = fv_tC/fv_dOF;
                if(fv_nC>5) continue;
                if(fv_nC < min_fvnC_2mu1tk){

                    iTrk=itk;
                    n_reco=2;

                    if(m_1.p()>m_2.p()){
                      j1=idxrec[i]; j2=idxrec[j];
                    }
                    else {j1=idxrec[j]; j2=idxrec[i];}
                    min_fvnC_2mu1tk = fv_nC;
                }

                if(!3muFitVtx.isValid()) continue;

                double fvnC2_tmp = fv.totalChiSquared()/fv.degreesOfFreedom();

                if (fvnC2_tmp > 2muTrkFitVtx_nC2_cut ) continue; // Store all the good candidates

                // Fill 2 muon info

                for (size_t it=0; it<2; ++it){

                    tmp_2mutrk_cand.mu_p[it] = sortedMu[it].p();
                    tmp_2mutrk_cand.mu_pt[it] = sotredMu[it].pt();
                    tmp_2mutrk_cand.mu_pterr[it] = sortedMu[it].muonBestTrack()->ptError();
                    tmp_2mutrk_cand.mu_bestPt[it] = sortedMu[it].muonBestTrack->pt();
                    tmp_2mutrk_cand.mu_eta[it] = sotredMu[it].eta();
                    tmp_2mutrk_cand.mu_phi[it] = sotredMu[it].phi();
                    tmp_2mutrk_cand.mu_charge[it] = sortedMu[it].charge();
                    tmp_2mutrk_cand.mu_vx[it] = sortedMu[it].vx();
                    tmp_2mutrk_cand.mu_vy[it] = sortedMu[it].vy();
                    tmp_2mutrk_cand.mu_vz[it] = sortedMu[it].vz();

                    // muon ID
                    bool _isGlobal = sortedMu[it].isGlobalMuon();
                    bool _isTracker = sortedMu[it].isTrackerMuon();
                    bool _isTrackerArb = muon::isGoodMuon(sortedMu[it].TrackerArbitrated);
                    bool _isRPC = sortedMu[it].isRPCMuon();
                    bool _isStandAlone = sortedMu[it].isStandAloneMuon();
                    bool _isPF = sortedMu[it].isPFMuon();

                    bool _hasInnerTrack = !mu.innerTrack().isNull();
                    bool _hasTunePTrack = !mu.tunePMuonBestTrack().isNull();
                    bool _hasPickyTrack = !mu.pickyTrack().isNull();
                    bool _hasDytTrack = !mu.dytTrack().isNull();
                    bool _hasTpfmsTrack = !mu.tpfmsTrack().isNull();

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
                    tmp_2mutrk_cand.mu_tRC2[it] = sortedMu[it].combinedQuality().trkRelChi2;
                    tmp_2mutrk_cand.mu_sRC2[it] = sortedMu[it].combinedQuality().staRelChi2;
                    tmp_2mutrk_cand.mu_Chi2LP[it] = sortedMu[it].combinedQuality().chi2LocalPosition;
                    tmp_2mutrk_cand.mu_Chi2LM[it] = sortedMu[it].combinedQuality().chi2LocalMomentum;
                    tmp_2mutrk_cand.mu_localDist[it] = sortedMu[it].combinedQuality().localDistance;
                    tmp_2mutrk_cand.mu_glbDEP[it] = sortedMu[it].combinedQuality().globalDeltaEtaPhi;
                    tmp_2mutrk_cand.mu_tightMatch[it] = sortedMu[it].combinedQuality().tightMatch;
                    tmp_2mutrk_cand.mu_glbTrkProb[it] = sortedMu[it].combinedQuality.glbTrackProbability;
                    tmp_2mutrk_cand.mu_calEM[it] = sortedMu[it].calEnergy().em;
                    tmp_2mutrk_cand.mu_calEMS9[it] = sortedMu[it].calEnergy().emS9;
                    tmp_2mutrk_cand.mu_calEMS25[it] = sortedMu[it].calEnergy().emS25;
                    tmp_2mutrk_cand.mu_calHad[it] = sortedMu[it].calEnergy().had;
                    tmp_2mutrk_cand.mu_calHadS9[it] = sortedMu[it].calEnergy.hadS9;;
                    tmp_2mutrk_cand.mu_nOMS[it] = sortedMu[it].numberOfMatchedStations();
                    tmp_2mutrk_cand.mu_nOM[it] = sortedMu[it].numberOfMatches(reco::Muon::segmentArbitraction);
                    tmp_2mutrk_cand.mu_comp2d[it] = muon::isGoodMuon(sortedMu[it], muon::TM2DCompatibilityTight);
                    tmp_2mutrk_cand.mu_calocomp[it] = muon::caloCompatibility(sortedMu[it]);
                    tmp_2mutrk_cand.mu_segcomp[it] = muon::segmentCompatibility(sortedMu[it]);

                    // Muon outer track info
                    tmp_2mutrk_cand.p_out[it] = _isStandAlone ? sortedMu[it].outerTrack()->p():-999;
                    tmp_2mutrk_cand.eta_out[it] = _isStandAlone ? sortedMu[it].outerTrack()->eta():-999;
                    tmp_2mutrk_cand.phi_out[it] = _isStandAlone ? sortedMu[it].outerTrack()->phi(): -999;
                    tmp_2mutrk_cand.outerchi2_reco[it] = _isStandAlone ? sortedMu[it].outerTrack()->normalizedChi2():-999;
                    tmp_2mutrk_cand.mu_glb.MSWVH_[it] = _isStandAlone ? sortedMu[it].outerTrack()->hitPattern().muonStationsWithValidHits():-999;
                    tmp_2mutrk_cand.qprod_reco[it] = _isStandAlone ? sortedMu[it].outerTrack()->charge()*sortedMu[it].innerTrack()->charge():-999;

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
                    tmp_2mutrk_cand.mu_inTrk_trkLayWithMeas[it] = _hasInnerTrack ? sortedMu[it].innerTracker()->hitPattern().trackerLayersWithMeasurement():-999;
                    tmp_2mutrk_cand.mu_inTrk_pixLayWithMeas[it] = _hasInnerTrack ? sortedMu[it].innerTracker()->hitPattern().pixelLayersWithMeasurement():-999;
                    tmp_2mutrk_cand.mu_inTrk_nValidTrkHits[it] = _hasInnerTrack ? sortedMu[it].innerTracker()->hitPattern().numberOfValidTrackerHits():-999;
                    tmp_2mutrk_cand.mu_inTrk_nValidPixHits[it] = _hasInnerTrack ? sortedMu[it].innerTracker()->hitPattern().numberOfValidTrackerHits():-999;
                    tmp_2mutrk_cand.mu_inTrk_validFraction[it] = _hasInnerTrack ? sortedMu[it].innerTracker()->validFraction():-999;
                    tmp_2mutrk_cand.mu_inTrk_nLostTrkHits[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfLostTrackerHits(hitPattern::TRACK_HITS):-999;
                    tmp_2mutrk_cand.mu_inTrk_nLostTrkHits_in[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfLostTrackerHits(hitPattern::MISSING_INNER_HITS):-999;
                    tmp_2mutrk_cand.mu_inTrk_nLostTrkHits_out[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->hitPattern().numberOfLostTrackerHits(hitPattern::MISSING_OUTER_HITS):-999;
                    tmp_2mutrk_cand.mu_inTrk_HP[it] = _hasInnerTrack ? sortedMu[it].innerTrack()->quality(TrackBase::highPurity);

                    // Isolation variables
                    reco::MuonIsolation IsoR03 = sortedMu[it].isolationR03;

                    tmp_2mutrk_cand.mu_IsoR03_sumPt[it] = IsoR03.sumPt;
                    tmp_2mutrk_cand.mu_IsoR03_nTrks[it] = IsoR03.nTracks;
                    tmp_2mutrk_cand.mu_IsoR03_emEt[it] = IsoR03.emEt;
                    tmp_2mutrk_cand.mu_IsoR03_hadEt[it] = IsoR03.hadEt;
                    tmp_2mutrk_cand.mu_IsoR03_emVetoEt[it] = IsoR03.emVetoEt;
                    tmp_2mutrk_cand.mu_IsoR03_hadVetoEt[it] = IsoR03.hadVetoEt;

                    // PF isolation variables
                    reco::MuonPFIsolation pfIso04 = sortedMu[it].pfIsolationR04();
                    reco::MuonPFIsolation pfIso03 = sortedMu[it].pfIsolationR03();

                    tmp_2mutrk_cand.chargedHadronIso   = pfIso04.sumChargedHadronPt;
                    tmp_2mutrk_cand.chargedHadronIsoPU = pfIso04.sumPUPt; 
                    tmp_2mutrk_cand.neutralHadronIso   = pfIso04.sumNeutralHadronEt;
                    tmp_2mutrk_cand.photonIso          = pfIso04.sumPhotonEt;

                    tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(sortedMu[it].pt(),sortedMu[it].eta(),sortedMu[it].phi(),
                            sortedMu[it].charge(),sortedMu[it].muonBestTrack()->ptError()));

                    tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(hasInnerTrack ? sortedMu[it].innerTrack()->pt()  : -1000.,
                            hasInnerTrack ? sortedMu[it].innerTrack()->eta() : -1000.,
                            hasInnerTrack ? sortedMu[it].innerTrack()->phi() : -1000.,
                            hasInnerTrack ? sortedMu[it].innerTrack()->charge()  : -1000.,
                            hasInnerTrack ? sortedMu[it].innerTrack()->ptError() : -1000.));

                    tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(isStandAlone ? sortedMu[it].outerTrack()->pt()  : -1000.,
                            isStandAlone ? sortedMu[it].outerTrack()->eta() : -1000.,
                            isStandAlone ? sortedMu[it].outerTrack()->phi() : -1000.,
                            isStandAlone ? sortedMu[it].outerTrack()->charge()  : -1000.,
                            isStandAlone ? sortedMu[it].outerTrack()->ptError() : -1000.));

                    tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(isGlobal ? sortedMu[it].globalTrack()->pt()  : -1000.,
                            isGlobal ? sortedMu[it].globalTrack()->eta() : -1000.,
                            isGlobal ? sortedMu[it].globalTrack()->phi() : -1000.,
                            isGlobal ? sortedMu[it].globalTrack()->charge()  : -1000.,
                            isGlobal ? sortedMu[it].globalTrack()->ptError() : -1000.));

                    tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->pt()  : -1000.,
                            hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->eta() : -1000.,
                            hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->phi() : -1000.,
                            hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->charge()  : -1000.,
                            hasTunePTrack ? sortedMu[it].tunePMuonBestTrack()->ptError() : -1000.));

                    tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(hasPickyTrack ? sortedMu[it].pickyTrack()->pt()  : -1000.,
                            hasPickyTrack ? sortedMu[it].pickyTrack()->eta() : -1000.,
                            hasPickyTrack ? sortedMu[it].pickyTrack()->phi() : -1000.,
                            hasPickyTrack ? sortedMu[it].pickyTrack()->charge()  : -1000.,
                            hasPickyTrack ? sortedMu[it].pickyTrack()->ptError() : -1000.));

                    tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(hasDytTrack ? sortedMu[it].dytTrack()->pt()  : -1000.,
                            hasDytTrack ? sortedMu[it].dytTrack()->eta() : -1000.,
                            hasDytTrack ? sortedMu[it].dytTrack()->phi() : -1000.,
                            hasDytTrack ? sortedMu[it].dytTrack()->charge()  : -1000.,
                            hasDytTrack ? sortedMu[it].dytTrack()->ptError() : -1000.));

                    tmp_2mutrk_cand.fits.push_back(tau23mu::MuonFit(hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->pt()  : -1000.,
                            hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->eta() : -1000.,
                            hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->phi() : -1000.,
                            hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->charge()  : -1000.,
                            hasTpfmsTrack ? sortedMu[it].tpfmsTrack()->ptError() : -1000.));

                    tmp_2mutrk_cand.mu_dxyBest  = -999; 
                    tmp_2mutrk_cand.mu_dzBest   = -999; 
                    tmp_2mutrk_cand.mu_dxyInner = -999; 
                    tmp_2mutrk_cand.mu_dzInner  = -999; 

                    tmp_2mutrk_cand.mu_dxybs = _isGlobal ? mu.globalTrack()->dxy(beamSpot->position()) :
                      _hasInnerTrack ? mu.innerTrack()->dxy(beamSpot->position()) : -1000;
                    tmp_2mutrk_cand.mu_dzbs  = _isGlobal ? mu.globalTrack()->dz(beamSpot->position()) :
                      _hasInnerTrack ? mu.innerTrack()->dz(beamSpot->position()) : -1000;

                    double dxy = -1000.;
                    double dz  = -1000.;

                    tmp_2mutrk_cand.isSoft    = 0;   
                    tmp_2mutrk_cand.isTight   = 0;   
                    tmp_2mutrk_cand.isHighPt  = 0;
                    tmp_2mutrk_cand.isLoose   = muon::isLooseMuon(sortedMu[it])  ? 1 : 0;   
                    tmp_2mutrk_cand.isMedium  = muon::isMediumMuon(sortedMu[it]) ? 1 : 0;   
                    if ( mu.isMatchesValid() && _isTrackerArb )
                    {
                      for ( reco::MuonChamberMatch match : mu.matches() )
                      {
                        muon_pog::ChambMatch ntupleMatch;

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

                            tmp_2mutrk_cand.matches.push_back(ntupleMatch);
                        }
                      }
                    }
                    tmp_2mutrk_cand.mu_isoPflow04 = (pfIso04.sumChargedHadronPt+ 
                        std::max(0.,pfIso04.sumPhotonEt+pfIso04.sumNeutralHadronEt - 0.5*pfIso04.sumPUPt)) / mu.pt();

                    tmp_2mutrk_cand.mu_isoPflow03 = (pfIso03.sumChargedHadronPt+ 
                        std::max(0.,pfIso03.sumPhotonEt+pfIso03.sumNeutralHadronEt - 0.5*pfIso03.sumPUPt)) / mu.pt();

                    if (vertexes->size() > 0) {
                      const reco::Vertex & vertex = vertexes->at(0);

                      tmp_2mutrk_cand.mu_dxy[it] = _isGlobal ? sortedMu[it].globalTrack()->dxy(vertex.position()) :
                        _hasInnerTrack ? sortedMu[it].innerTrack()->dxy(vertex.position()) : -1000;
                      tmp_2mutrk_cand.mu_dz[it] = _isGlobal ? sortedMu[it].globalTrack()->dz(vertex.position()) :
                        _hasInnerTrack ? sortedMu[it].innerTrack()->dz(vertex.position()) : -1000;

                      tmp_2mutrk_cand.mu_dxyBest[it]  = sortedMu[it].muonBestTrack()->dxy(vertex.position()); 
                      tmp_2mutrk_cand.mu_dzBest[it]   = sortedMu[it].muonBestTrack()->dz(vertex.position()); 
                      if(_hasInnerTrack) { 
                        tmp_2mutrk_cand.mu_dxyInner[it] = sortedMu[it].innerTrack()->dxy(vertex.position()); 
                        tmp_2mutrk_cand.mi_dzInner[it]  = sortedMu[it].innerTrack()->dz(vertex.position()); 
                      } 

                      tmp_2mutrk_cand.isSoft[it]    = muon::isSoftMuon(sortedMu[it],vertex)   ? 1 : 0;   
                      tmp_2mutrk_cand.isTight[it]   = muon::isTightMuon(sortedMu[it],vertex)  ? 1 : 0;   
                      tmp_2mutrk_cand.isHighPt[it]  = muon::isHighPtMuon(sortedMu[it],vertex) ? 1 : 0;

                    }
                    if(sortedMu[it].isTimeValid()) { 
                      tmp_2mutrk_cand.mu_TimeDof = sortedMu[it].time().nDof; 
                      tmp_2mutrk_cand.mu_Time    = sortedMu[it].time().timeAtIpInOut; 
                      tmp_2mutrk_cand.mu_TimeErr = sortedMu[it].time().timeAtIpInOutErr; 
                    } 
                    else { 
                      tmp_2mutrk_cand.muonTimeDof = -999; 
                      tmp_2mutrk_cand.muonTime    = -999; 
                      tmp_2mutrk_cand.muonTimeErr = -999; 
                    } 

                    if(sortedMu[it].rpcTime().nDof > 0) { 
                      tmp_2mutrk_cand.mu_RpcTimeDof[it] = sortedMu[it].rpcTime().nDof; 
                      tmp_2mutrk_cand.mu_RpcTime[it]    = sortedMu[it].rpcTime().timeAtIpInOut; 
                      tmp_2mutrk_cand.mu_RpcTimeErr[it] = sortedMu[it].rpcTime().timeAtIpInOutErr; 
                    } 
                    else { 
                      tmp_2mutrk_cand.mu_RpcTimeDof[it] = -999; 
                      tmp_2mutrk_cand.mu_RpcTime[it]   = -999; 
                      tmp_2mutrk_cand.mu_RpcTimeErr[it] = -999; 
                    } 


                    tmp_2mutrk_cand.mu_dxy[it]  = _isGlobal ? sortedMu[it].globalTrack()->dxy(vertex.position());
                    tmp_2mutrk_cand.mu_dz[it]  = _isGlobal ? sortedMu[it].globalTrack()->dxy(vertex.position());
                    tmp_2mutrk_cand.mu_edxy[it]   = _isGlobal ? sortedMu[it].globalTrack()->dxyError() : _hasInnerTrack ? sortedMu[it].innerTrack()->dxyError() : -1000;
                    tmp_2mutrk_cand.mu_edz[it]    = _isGlobal ? sortedMu[it].globalTrack()->dzError()  : _hasInnerTrack ? sortedMu[it].innerTrack()->dzError() : -1000;

                    tmp_2mutrk_cand.mu_TMOST[it] = muon::isGoodMuon(mu[i], muon::TMOneStationTight);
                    tmp_2mutrk_cand.mu_TMOSAT[it] = muon::isGoodMuon(mu[i], muon::TMOneStationAngTight);
                    tmp_2mutrk_cand.mu_TMLST[it] = muon::isGoodMuon(mu[i], muon::TMLastStationTight);
                    tmp_2mutrk_cand.mu_TMLSAT[it] = muon::isGoodMuon(mu[i], muon::TMLastStationAngTight);
                    tmp_2mutrk_cand.mu_TMLSOLPT[it] = muon::isGoodMuon(mu[i], muon::TMLastStationOptimizedLowPtTight);
                    tmp_2mutrk_cand.mu_TMLSOBLPT[it] = muon::isGoodMuon(mu[i], muon::TMLastStationOptimizedBarrelLowPtTight);

                    // Tau-Muon variables
                    mu_phi_dR[it] = deltaR(sortedMu[it].eta(), sortedMu[it].phi(), vtau.Eta(), vtau.Phi());
                }

                // Fill track info

                tmp_2mutrk_cand.pterr = trk->ptError();
                tmp_2mutrk_cand.trkhp = trk->quality(TrackBase::highPurity);
                tmp_2mutrk_cand.charge = trk->charge();
                tmp_2mutrk_cand.trk_p = trk->p();
                tmp_2mutrk_cand.trk_pt = trk->pt();
                tmp_2mutrk_cand.trk_eta = trk->eta();
                tmp_2mutrk_cand.trk_phi = trk->phi();
                tmp_2mutrk_cand.nOVPH = trk->hitPattern().numberOfValidPixelHits();
                tmp_2mutrk_cand.iTvF = trk->validFraction();
                tmp_2mutrk_cand.tLWM = trk->hitPattern().trackerLayersWithMeasurement();
                tmp_2mutrk_cand.pLWM = trk->hitPattern().pixelLayersWithMeasurement();
                tmp_2mutrk_cand.nOVTH = trk->hitPattern().numberOfValidTrackerHits();
                tmp_2mutrk_cand.nOLTH = trk->hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS);
                tmp_2mutrk_cand.nOLTHin = trk->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);
                tmp_2mutrk_cand.nOLTHout = trk->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS);
                tmp_2mutrk_cand.iTnC = trk->normalizedChi2();

                // Compute 2mu + trk invaiant mass
                TLorentzVector vec_mu1, vec_mu2, vec_trk;

                vec_mu1.SetPtEtaPhiM(sortedMu[0].pt(), sortedMu[0].eta(), sortedMu[0].phi(), MU_MASS);
                vec_mu2.SetPtEtaPhiM(sortedMu[1].pt(), sortedMu[1].eta(), sortedMu[1].phi(), MU_MASS);
                vec_trk.SetPtEtaPhiM(trk.pt(), trk.eta(), trk.phi(), PI_MASS);

                2muTrk_InvMass.push_back((vec_mu1+vec_mu2+vec_mu3).M());
                2muTrkVtx_nC2.push_back(fvnC2_tmp);
					 
					 // Refit tracks 
                
					 TrackRef trk1 = sortedMu[0].innerTrack();
                TrackRef trk2 = sortedMu[1].innerTrack();
                ESHandle<TransientTrackBuilder> theB;
                iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

                // calculate closest approach ===> format has to change make a class and then get all the output results
                ClosestApproachInRPhi cApp12, cApp23, cApp31;
                cApp12.calculate(t_trks[0].initialFreeState(), t_trks[1].initialFreeState());
                cApp23.calculate(t_trks[1].initialFreeState(), t_trks[2].initialFreeState());
                cApp31.calculate(t_trks[2].initialFreeState(), t_trks[0].initialFreeState());

                if(!(cApp12.status()&&cApp23.status()&&cApp31.status())) { return; cout<<"DCA unvalid!"<<endl; }
                tmp_2mutrk_cand.dca12 = cApp12.distance();
                tmp_2mutrk_cand.dca23_reco = cApp23.distance();
                tmp_2mutrk_cand.dca31_reco = cApp31.distance();
                tmp_2mutrk_cand.dca_max = TMath::Max(dca12_reco, TMath::Max(dca23_reco, dca31_reco));

                TransientVertex fv = tau23mu::fitTriMuVertex(trk1, trk2, trk);
                if(!fv.isValid()) { return; cout<<"Vertex Fit unvalid!"<<endl; }

                TLorentzVector vtau_refit, vmu_refit;
                vtau_refit.SetPtEtaPhiM(0, 0, 0, 0);
                vector<TransientTrack>::const_iterator trkIt = fv.refittedTracks().begin();
                for(; trkIt != fv.refittedTracks().end(); ++ trkIt) {
                    const Track & trkrefit = trkIt->track();
                    vmu_refit.SetPtEtaPhiM(trkrefit.pt(), trkrefit.eta(), trkrefit.phi(), (n_reco>2?0.106:0.140));
                    vtau_refit += vmu_refit;
                }

                m3mu_refit = vtau_refit.M();

                fv_tC = fv.totalChiSquared();
                int fv_dOF = fv.degreesOfFreedom();
                fv_nC = fv_tC/fv_dOF;
                fv_Prob = TMath::Prob(fv_tC,(int)fv_dOF);

                vector<TransientTrack> t_trks12, t_trks23, t_trks31;
                t_trks12.push_back(theB->build(trk1)); t_trks12.push_back(theB->build(trk2));
                t_trks23.push_back(theB->build(trk2)); t_trks23.push_back(theB->build(t3));
                t_trks31.push_back(theB->build(t3)); t_trks31.push_back(theB->build(trk1));
                KalmanVertexFitter kvf_trks12, kvf_trks23, kvf_trks31;
                TransientVertex fv_trks12 = kvf_trks12.vertex(t_trks12);
                TransientVertex fv_trks23 = kvf_trks23.vertex(t_trks23);
                TransientVertex fv_trks31 = kvf_trks31.vertex(t_trks31);

                tmp_2mutrk_cand.fvwo_tC[0] = fv_trks23.totalChiSquared();
                tmp_2mutrk_cand.fvwo_nC[0] = fvwo_tC[0]/fv_trks23.degreesOfFreedom();
                tmp_2mutrk_cand.fvwo_tC[1] = fv_trks31.totalChiSquared();
                tmp_2mutrk_cand.fvwo_nC[1] = fvwo_tC[1]/fv_trks31.degreesOfFreedom();
                tmp_2mutrk_cand.fvwo_tC[2] = fv_trks12.totalChiSquared();
                tmp_2mutrk_cand.fvwo_nC[2] = fvwo_tC[2]/fv_trks12.degreesOfFreedom();
              }
            }
        }

        TVector3 vtauxyz(vtau.Px(), vtau.Py(), vtau.Pz());

        ////////////////////
        // find the good PV
        double PVZ = fv.position().z()-fv_dxy*vtau.Pz()/vtau.Pt();
        double dispv1 = 99, dispvgen=99, dphi_pv = -1;
        //int ipvPVZ = 99, ipvgen = 99;
        for(size_t jpv = 0; jpv < pvs->size(); jpv++) {
            const Vertex & vi = (*pvs)[jpv];

            if(abs(vi.position().Z()-PVZ)<dispv1){
              dispv1=abs(vi.position().Z()-PVZ);
              ipv1=jpv;
            }
            if(abs(vi.position().Z()-gen_pv)<dispvgen){
              dispvgen=abs(vi.position().Z()-gen_pv);
              ipv_gen=jpv;
            }
            TVector3 Dv3D_reco(fv.position().x() - vi.x(), fv.position().y() - vi.y(), fv.position().z() - vi.z());
            double Cosdphi_3D = Dv3D_reco.Dot(vtauxyz)/(Dv3D_reco.Mag()*vtauxyz.Mag());
            if(Cosdphi_3D>dphi_pv){
              dphi_pv = Cosdphi_3D;
              ipv2=jpv;
            }

        }
        const Vertex & pv0 = (*pvs)[ipv2];

        //////////////////////////////////////////////////
        // refit PV with and w.o. the 3 mu
        //////////////////////////////////////////////////
        
		  pv1_tC = 999; pv1_nC = 99; pv2_tC = 999; pv2_nC = 99;
        vector<TransientTrack> pv_trks;
        TransientVertex pv2, pv1;

        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
        for(Vertex::trackRef_iterator itk = pv0.tracks_begin(); itk != pv0.tracks_end(); itk++) {
            if(trk.pt()>1) {
              if(deltaR(sortedMu[0].eta(), sortedMu[0].phi(), trk.eta(), trk.phi())<0.01)continue;
              if(deltaR(sortedMu[1].eta(), sortedMu[1].phi(), trk.eta(), trk.phi())<0.01)continue;
              if(deltaR(trk->eta(), trk->phi(), trk.eta(), trk.phi())<0.01)continue;
            }
            pv_trks.push_back(theB->build(trk));
        }
        if(pv_trks.size()>1) {
            KalmanVertexFitter kvf_pv;
            pv1 = kvf_pv.vertex(pv_trks);
            if(pv1.isValid()){
              pv1_tC = pv1.totalChiSquared();
              pv1_nC = pv1_tC/pv1.degreesOfFreedom();
            }

            // adding the 3 mu tracks
            pv_trks.push_back(theB->build(trk1));
            pv_trks.push_back(theB->build(trk2));
            pv_trks.push_back(theB->build(t3));
            pv2 = kvf_pv.vertex(pv_trks);
            if(pv2.isValid()){
              pv2_tC = pv2.totalChiSquared();
              pv2_nC = pv2_tC/pv2.degreesOfFreedom();
            }
        }

        Vertex pvv = pv0;  // the final PV
        if(pv1.isValid()) pvv = Vertex(pv1);
        math::XYZPoint pv1P = math::XYZPoint(pvv.x(), pvv.y(), pvv.z());

        d0_reco[0] = abs(mu[0].innerTrack()->dxy(pv1P));
        d0_reco[1] = abs(mu[1].innerTrack()->dxy(pv1P));
        d0_reco[2] = abs(t3->dxy(pv1P));
        d0sig_reco[0] = -1; d0sig_reco[1] = -1; d0sig_reco[2] = -1;
        GlobalVector dir1(mu[0].px(), mu[0].py(), mu[0].pz());
        GlobalVector dir2(mu[1].px(), mu[1].py(), mu[1].pz());
        GlobalVector dir3(t3->px(), t3->py(), t3->pz());
        std::pair<bool, Measurement1D> ip2d_1 = IPTools::signedTransverseImpactParameter(t_trks[0], dir1, pvv);
        std::pair<bool, Measurement1D> ip2d_2 = IPTools::signedTransverseImpactParameter(t_trks[1], dir2, pvv);
        std::pair<bool, Measurement1D> ip2d_3 = IPTools::signedTransverseImpactParameter(t_trks[2], dir3, pvv);
        if(ip2d_1.first) d0sig_reco[0] = abs(ip2d_1.second.value()/ip2d_1.second.error());
        if(ip2d_2.first) d0sig_reco[1] = abs(ip2d_2.second.value()/ip2d_2.second.error());
        if(ip2d_3.first) d0sig_reco[2] = abs(ip2d_3.second.value()/ip2d_3.second.error());

        ////////////////////
        // displacement 2D
        TVector3 dv2D_reco(fv.position().x() - pv1P.x(), fv.position().y() - pv1P.y(), 0);
        TVector3 vtauxy(vtau.Px(), vtau.Py(), 0);
        fv_cosdphi = dv2D_reco.Dot(vtauxy)/(dv2D_reco.Perp()*vtauxy.Perp());
        VertexDistanceXY vdistXY;
        Measurement1D distXY = vdistXY.distance(Vertex(fv), pvv);
        fv_dxy = distXY.value();
        fv_dxysig = distXY.significance();
        fv_ppdl = distXY.value()*fv_cosdphi * m3mu_reco/vtauxy.Perp();

        ////////////////////
        // displacement 3D
        TVector3 dv3D_reco(fv.position().x() - pv1P.x(), fv.position().y() - pv1P.y(), fv.position().z() - pv1P.z());
        //TVector3 vtauxyz(vtau.Px(), vtau.Py(), vtau.Pz());
        fv_cosdphi3D = dv3D_reco.Dot(vtauxyz)/(dv3D_reco.Mag()*vtauxyz.Mag());
        VertexDistance3D dist;
        fv_d3D = dist.distance(Vertex(fv), pvv).value(); // = dv_reco.Mag() ??
        fv_d3Dsig = dist.distance(Vertex(fv), pvv).significance();
        fv_ppdl3D = fv_d3D*fv_cosdphi3D*m3mu_reco/vtau.P();
        vector<double> softmueta, softmuphi;
        pv_nmu = 0;

        for(size_t i = 0; i < muons->size(); ++ i) {
            if(i==j1)continue;
            if(i==j2)continue;
            if(i==j3)continue;
            const Muon & m_1 = (*muons)[i];
            if(!(abs(m_1.eta())<2.4)) continue;
            if(!(muon::isGoodMuon(m_1, muon::TMOneStationTight))) continue;
            if(!(m_1.innerTrack()->hitPattern().trackerLayersWithMeasurement()>5))continue;
            if(!(m_1.innerTrack()->hitPattern().pixelLayersWithMeasurement()>0))continue;
            //if(!(abs(m_1.innerTrack()->dxy(pv0.position())) < .3))continue;
            if(!(abs(m_1.innerTrack()->dz(pv1P)) < 1))continue;
            pv_nmu++;
            softmueta.push_back(m_1.eta());
            softmuphi.push_back(m_1.phi());
        }

        ////////////////////
        // secondary vertices
        n_sv = 0 ;
        for(size_t isv = 0; isv < svs->size(); isv++) {
            const Vertex & sv = (*svs)[isv];
            if(abs(sv.p4().M()-0.498)<.03 && sv.tracksSize()==2)continue; // no Ks

            double dx = sv.x()-pv1P.x();
            double dy = sv.y()-pv1P.y();
            double dz = sv.z()-pv1P.z();
            TVector3 sv_reco(dx, dy, dz);
            sv_overlap[n_sv]=deltaR(sv_reco.Eta(), sv_reco.Phi(), dv3D_reco.Eta(), dv3D_reco.Phi());

            TVector3 svxyz(sv.p4().Px(), sv.p4().Py(), sv.p4().Pz());
            sv_cosdphi3D[n_sv] = sv_reco.Dot(svxyz)/(sv_reco.Mag()*svxyz.Mag());
            VertexDistance3D distsv;
            sv_d3D[n_sv] = distsv.distance(sv, pvv).value();
            sv_d3Dsig[n_sv] = distsv.distance(sv, pvv).significance();
            sv_ppdl3D[n_sv] = sv_d3D[n_sv]*sv_cosdphi3D[n_sv]*sv.p4().M()/sv.p4().P();

            sv_nmu[n_sv] = 0;
            for(Vertex::trackRef_iterator itk = sv.tracks_begin(); itk != sv.tracks_end(); itk++) {
              for(size_t imu = 0; imu < softmueta.size(); imu ++) {
                if(deltaR(softmueta[imu], softmuphi[imu], (**itk).eta(), (**itk).phi())<0.01)
                    sv_nmu[n_sv] ++;
              }
            }

            sv_mass[n_sv] = sv.p4().M();
            sv_pt[n_sv] = sv.p4().Pt();
            sv_dz[n_sv] = abs(dz);
            sv_ntrk[n_sv] = sv.tracksSize();
            n_sv++;
        }

        dzpv_reco[0] = 1;
        dzpv_reco[1] = 1;
        dzpv_reco[2] = 1;

        dz_reco[0] = abs(mu[0].innerTrack()->dz(pv1P));
        dz_reco[1] = abs(mu[1].innerTrack()->dz(pv1P));
        dz_reco[2] = abs(t3->dz(pv1P));

        for(size_t jpv = 0; jpv < pvs->size(); jpv++) {
            if(jpv==ipv2)continue;
            const Vertex & vi = (*pvs)[jpv];
            if(abs(mu[0].innerTrack()->dz(vi.position()))<dz_reco[0]) dzpv_reco[0]=-1;
            if(abs(mu[1].innerTrack()->dz(vi.position()))<dz_reco[1]) dzpv_reco[1]=-1;
            if(abs(t3->dz(vi.position()))<dz_reco[2]) dzpv_reco[2]=-1;
        }

        ////////////////////
        // Track Isolation
        // How to decide if a track is associated with a certain PV ?
        double pttrk_tau = 0, pttrk_tau05 = 0,  pttrk_m1 = 0, pttrk_m2 = 0, pttrk_m3 = 0;
        mindca_iso = 99; mindca_iso05 = 99;
        ntrk_tau = 0; ntrk_tau05 = 0; ntrk_tau_b = 0; ntrk_sum = 0;
        ntrk_reco[0] = 0;  ntrk_reco[1] = 0;  ntrk_reco[2] = 0;
        ntrk0p1 = 0; ntrk0p2 = 0; ntrk0p5 = 0; maxdxy_pv0 =0;

        math::XYZPoint fvP = math::XYZPoint(fv.position().x(), fv.position().y(), fv.position().z());
        for(size_t i = 0; i < trks->size(); i++) {
            const Track & t = (*trks)[i];
            if(!(t.quality(TrackBase::tight)))continue;
            if(deltaR(mu[0].eta(), mu[0].phi(), t.eta(), t.phi())<0.01)continue;
            if(deltaR(mu[1].eta(), mu[1].phi(), t.eta(), t.phi())<0.01)continue;
            if(deltaR(t3->eta(), t3->phi(), t.eta(), t.phi())<0.01)continue;

            double dz = abs(t.dz(fvP));
            double dxy = abs(t.dxy(fvP));
            double dca_fv = sqrt(dz*dz+dxy*dxy);
            double dr_tau = deltaR(t.eta(), t.phi(), vtau.Eta(), vtau.Phi());

            // iso no. 1b - using pt_min, drtau_max of the 3 mu
            if(t.pt() > 0.33*pt_min && dr_tau < 3.*drtau_max && dca_fv<0.05 ) {
              pttrk_tau += t.pt();
              ntrk_tau++; // iso 3b
              if(dca_fv<mindca_iso)mindca_iso=dca_fv; // iso 4b
            } 

            if(t.pt()<1.0) continue;  // was 1.2
            // iso no. 1
            if(dr_tau < 0.5 && dca_fv<0.05 ) {
              pttrk_tau05 += t.pt();
              ntrk_tau05++; // iso 3
              //if(dca_fv<mindca_iso05)mindca_iso05=dca_fv; // iso 4
            }

            if(dca_fv<0.05)ntrk_tau_b++; // iso 3b
            if(dca_fv<mindca_iso05)mindca_iso05=dca_fv; // iso 4

            TransientTrack trkiso = theB->build(t);
            ClosestApproachInRPhi cAppm1, cAppm2, cAppm3;
            cAppm1.calculate(trkiso.initialFreeState(), t_trks[0].initialFreeState());
            cAppm2.calculate(trkiso.initialFreeState(), t_trks[1].initialFreeState());
            cAppm3.calculate(trkiso.initialFreeState(), t_trks[2].initialFreeState());
            if(!(cAppm1.status()&&cAppm2.status()&&cAppm3.status())) continue;

            // iso no. 2
            if(deltaR(t.eta(), t.phi(), mu[0].eta(), mu[0].phi()) < 0.3 && cAppm1.distance() < 0.1) {// && dz1 < .3) 
              ntrk_reco[0]++;
              pttrk_m1 += t.pt();
            }
            if(deltaR(t.eta(), t.phi(), mu[1].eta(), mu[1].phi()) < 0.3 && cAppm2.distance() < 0.1) {//&& dz2 < .3) 
              ntrk_reco[1]++;
              pttrk_m2 += t.pt();
            }
            if(deltaR(t.eta(), t.phi(), t3->eta(), t3->phi()) < 0.3 && cAppm3.distance() < 0.1) {//&& dz3 < .3) 
              ntrk_reco[2]++;
              pttrk_m3 += t.pt();
            }
            if( (deltaR(t.eta(), t.phi(), mu[0].eta(), mu[0].phi()) < 0.3 && cAppm1.distance() < 0.1 )
                ||(deltaR(t.eta(), t.phi(), mu[1].eta(), mu[1].phi()) < 0.3 && cAppm2.distance() < 0.1 )
                ||(deltaR(t.eta(), t.phi(), t3->eta(), t3->phi()) < 0.3 && cAppm3.distance() < 0.1 )
            ) ntrk_sum++;


            // displaced track counting
            // only tracks consistent with PV
            double dz_pv0=abs(t.dz(pv1P));
            if(!(dz_pv0 < 1))continue;
            double dxy_pv0 = abs(t.dxy(pv1P));
            if(dxy_pv0>0.1) ntrk0p1++;
            if(dxy_pv0>0.2) ntrk0p2++;
            if(dxy_pv0>0.5) ntrk0p5++;
            if(dxy_pv0>maxdxy_pv0) maxdxy_pv0 = dxy_pv0;

        }

        trkrel_tau = pttrk_tau/vtau.Pt();
        trkrel_tau05 = pttrk_tau05/vtau.Pt();
        trkrel_reco[0] = pttrk_m1/mu[0].pt(); trkrel_reco[1] = pttrk_m2/mu[1].pt(); trkrel_reco[2] = pttrk_m3/t3->pt();
        trkrel_max = TMath::Max(trkrel_reco[0], TMath::Max( trkrel_reco[1], trkrel_reco[2]));


        // Good Global Muon, Tight Global Muon
        ggm_reco[0] = (glbnC_reco[0]<3 && cLP_reco[0]<12 && tKink_reco[0]<20 && segmcomp_reco[0]>0.303)?1:0;
        tgm_reco[0] = (glbnC_reco[0]<10 && pf_reco[0] && nOVMH_reco[0]>0 && nOMS_reco[0]>1
              && d0_reco[0]<0.2 && dz_reco[0]<0.5 && nOVPH_reco[0]>0 && tLWM_reco[0]>5)?1:0;

      }



    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Tau23MuNtupleMaker);
