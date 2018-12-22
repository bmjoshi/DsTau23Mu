#ifndef 3MU_CANDIDATE
#define 3MU_CANDIDATE

namespace tau3mu{
   public:

    class 3mu_cand(){
      public:
        virtual void reset(){

          for (size_t i=0; i<3; i++){
           mu_pt[i] = 0;
           mu_eta[i] = 0;
           mu_phi[i] = 0;
           mu_charge[i] = 0;
          }
          mu_pt_max = 0;
          mu_eta_max = 0;
          mu_pt_min = 0;
          mu_eta_min = 0;

        }

        Double_t mu_pt_max;
        Double_t mu_pt_min;
        Double_t mu_eta_max;
        Double_t mu_eta_min;
        
		  std::array<Double_t, 3> mu_p;
        std::array<Double_t, 3> mu_pt;
        std::array<Double_t, 3> mu_eta;
        std::array<Double_t, 3> mu_phi;
        std::array<Double_t, 3> mu_charge;
        std::array<Double_t, 3> mu_uSta;
        std::array<Double_t, 3> mu_trkKink;
        std::array<Double_t, 3> mu_log_glbKink;
        std::array<Double_t, 3> mu_trkRelChi2;
        std::array<Double_t, 3> mu_staRelChi2;
        std::array<Double_t, 3> mu_Chi2LP;
        std::array<Double_t, 3> mu_Chi2LM;
        std::array<Double_t, 3> mu_localDist;
        std::array<Double_t, 3> mu_staRelChi2;
        std::array<Double_t, 3> mu_glbDEP;
        std::array<Double_t, 3> mu_tightMatch;
        std::array<Double_t, 3> mu_glbTrkProb;
        std::array<Double_t, 3> calEM_reco;
        std::array<Double_t, 3> calEMS9_reco;
        std::array<Double_t, 3> calEMS25_reco;
        std::array<Double_t, 3> calHad_reco;
        std::array<Double_t, 3> calHadS9_reco;
        std::array<Int_t, 3> nOMS_reco;
        std::array<Int_t, 3> nOM_reco;
        std::array<Double_t, 3> comp2d_reco;
        std::array<Double_t, 3> calocomp_reco;
        std::array<Double_t, 3> trkhp_reco;
        std::array<Double_t, 3> segcomp_reco;
        std::array<Double_t, 3> trkHP_reco;
        
		  std::array<bool, 3> mu_isPF;
        std::array<bool, 3> mu_isRPC;
        std::array<bool, 3> mu_isStandAlone;
        std::array<bool, 3> mu_isTrackerArb;

		  std::array<bool, 3> mu_isTracker;
		  std::array<bool, 3> mu_isGlobal;
		  std::array<bool, 3> mu_isSoft;
		  std::array<bool, 3> mu_isTight;
		  std::array<bool, 3> mu_isLoose;
		  std::array<bool, 3> mu_isMedium;
		  std::array<bool, 3> mu_isHighPt;

        std::array<Double_t, 3> segcomp_reco;
        std::array<Double_t, 3> segcomp_reco;
        std::array<Double_t, 3> mu_ptErr;
        
		  std::array<Double_t, 3> mu_inTrk_p;
        std::array<Double_t, 3> mu_inTrk_pt;
        std::array<Double_t, 3> mu_inTrk_eta;
        std::array<Double_t, 3> mu_inTrk_phi;
        std::array<Int_t, 3> mu_inTrk_trkLayWithMeas;
        std::array<Int_t, 3> mu_inTrk_pixLayWithMeas;
        std::array<Int_t, 3> mu_inTrk_nValidTrkHits;
        std::array<Int_t, 3> mu_inTrk_nLostTrkHits;
        std::array<Int_t, 3> mu_inTrk_nLostTrkHits_in;
        std::array<Int_t, 3> mu_inTrk_nLostTrkHits_out;
		  
		  std::array<Double_t, 3> mu_tau_dR;
		  
		  std::array<Double_t, 3> mu_IsoR03_sumPt;
		  std::array<Double_t, 3> mu_IsoR03_nTrks;
		  std::array<Double_t, 3> mu_IsoR03_emEt;
		  std::array<Double_t, 3> mu_IsoR03_emEt;
		  std::array<Double_t, 3> mu_IsoR03_hadEt;
		  std::array<Double_t, 3> mu_IsoR03_emVetoEt;
		  std::array<Double_t, 3> mu_IsoR03_hadVetoEt;

        iTnC_reco[i] = mu[i].innerTrack()->normalizedChi2();
        tma_reco[i] = muon::isGoodMuon(mu[i], muon::TrackerMuonArbitrated);
        tmost_reco[i] = muon::isGoodMuon(mu[i], muon::TMOneStationTight);
        tmosat_reco[i] = muon::isGoodMuon(mu[i], muon::TMOneStationAngTight);
        tmlst_reco[i] = muon::isGoodMuon(mu[i], muon::TMLastStationTight);
        tmlsat_reco[i] = muon::isGoodMuon(mu[i], muon::TMLastStationAngTight);
        tmlsolpt_reco[i] = muon::isGoodMuon(mu[i], muon::TMLastStationOptimizedLowPtTight);
        tmlsoblpt_reco[i] = muon::isGoodMuon(mu[i], muon::TMLastStationOptimizedBarrelLowPtTight);
        timeatipinouterr_reco[i] = mu[i].time().timeAtIpInOutErr;

        Int_t  isGlobal;
        Int_t  isTracker;
        Int_t  isTrackerArb;
        Int_t  isRPC;
        Int_t  isStandAlone;
        Int_t  isPF;

        Int_t  isSoft;
        Int_t  isLoose;
        Int_t  isTight;
        Int_t  isMedium;
        Int_t  isHighPt;

        //Detector Based Isolation
        Float_t trackerIso;
        Float_t EMCalIso;
        Float_t HCalIso;

        // PF Isolation
        Float_t chargedHadronIso;
        Float_t chargedHadronIsoPU;
        Float_t photonIso;
        Float_t neutralHadronIso;


        Float_t isoPflow04; // PF isolation in dR<0.4 cone dBeta
        Float_t isoPflow03; // PF isolation in dR<0.3 cone dBeta

        Float_t dxy;     // signed transverse distance to primary vertex [cm]
        Float_t dz;      // signed longitudinal distance to primary vertex at min. transv. distance [cm]
        Float_t edxy;    // uncertainty on dxy [cm]
        Float_t edz;     // uncertainty on dz [cm]
        Float_t dxybs;    // signed transverse distance to beamspot [cm]
        Float_t dzbs;    // signed longitudinal distance to beamspot [cm]
        Float_t vx;
        Float_t vy;
        Float_t vz;

        Int_t  nHitsGlobal;
        Int_t  nHitsTracker;
        Int_t  nHitsStandAlone; 

        // Variables for ID 
        //  - General (Tight, HighPt, Soft) 
        Float_t glbNormChi2; 
        Float_t trkNormChi2; 
        Int_t  trkMuonMatchedStations; 
        Int_t  glbMuonValidHits; 
        Int_t  trkPixelValidHits; 
        Int_t  trkPixelLayersWithMeas; 
        Int_t  trkTrackerLayersWithMeas; 

        //  - HighPt 
        Float_t bestMuPtErr; 

        //  - Medium 
        Float_t trkValidHitFrac; 
        Float_t trkStaChi2; 
        Float_t trkKink; 
        Float_t muSegmComp; 

        //  - Soft 
        Int_t  isTrkMuOST; 
        Int_t  isTrkHP; 

        Float_t dxyBest; 
        Float_t dzBest; 
        Float_t dxyInner; 
        Float_t dzInner; 

        // Muon time 
        Float_t muonTimeDof; 
        Float_t muonTime; 
        Float_t muonTimeErr;

        // Muon time 
        Float_t muonRpcTimeDof; 
        Float_t muonRpcTime; 
        Float_t muonRpcTimeErr;

        std::vector<HitInfo> hits;
        std::vector<ChambMatch> matches;
        std::vector<MuonFit> fits;

        3mu_cand(){};
        virtual ~3mu_cand(){};

        inline Float_t fitPt( const tau3mu::MuonFitType type ) 
        {
          return fits.at(type).pt;
        };

        inline Float_t fitEta( const tau3mu::MuonFitType type ) 
        {
          return fits.at(type).eta;
        };

        inline Float_t fitPhi( const tau3mu::MuonFitType type ) 
        {
          return fits.at(type).phi;
        };

        inline Int_t fitCharge( const tau3mu::MuonFitType type ) 
        {
          return fits.at(type).charge;
        };

        inline Float_t fitPtErr( const tau3mu::MuonFitType type ) 
        {
          return fits.at(type).ptErr;
        };

        ClassDef(3mu_cand,4)


    };

#endif 
