#ifndef DsTau23Mu_Tau23MuNTau23MuTree_H
#define DsTau23Mu_Tau23MuNTau23MuTree_H

#include "TROOT.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <array>
#include <algorithm>

namespace tau23mu {
  
  enum MuonDetType { DT=0, CSC, RPC };

  class ChambMatch {
  public:
    Int_t r; // station/disk
    Int_t phi; // sector
    Int_t eta; // ring/wheel

    MuonDetType type;

    Float_t x; 
    Float_t y;

    Float_t dXdZ; 
    Float_t dYdZ;

    Float_t dx;  // 999999 if not matched with a segment (I think) 
    Float_t dy;  // 999999 if not matched with a segment (I think)

    Float_t dDxDz;  // 999999 if not matched with a segment (I think) 
    Float_t dDyDz;  // 999999 if not matched with a segment (I think)

    Float_t errxTk; 
    Float_t erryTk; 

    Float_t errDxDzTk; 
    Float_t errDyDzTk; 

    Float_t errxSeg;  // 999999 if not matched with a segment (I think)
    Float_t errySeg;  // 999999 if not matched with a segment (I think) 

    Float_t errDxDzSeg;  // 999999 if not matched with a segment (I think)
    Float_t errDyDzSeg;  // 999999 if not matched with a segment (I think) 

    ChambMatch(){};
    virtual ~ChambMatch(){};

    ClassDef(ChambMatch,2)
  };

  class HitInfo {
  public:
    Int_t r; // station/disk
    Int_t phi; // sector
    Int_t eta; // ring/wheel

    MuonDetType type;

    Int_t nHits; 
    Int_t nHitsPhi; 
    Int_t nHitsTheta; 

    HitInfo(){};
    virtual ~HitInfo(){};

    ClassDef(HitInfo,1)
  };

  enum MuonFitType { DEFAULT=0, INNER, STA, GLB, TUNEP, PICKY, DYT, TPFMS};

  class MuonFit {
  public:
    Float_t pt;  // pt [GeV] 
    Float_t eta; // eta
    Float_t phi; // phi

    Int_t charge;  // charge

    Float_t ptErr; // fit sigma pT 

    MuonFit(){};
    MuonFit(Float_t in_pt,
      Float_t in_eta,
      Float_t in_phi,
      Int_t in_charge,
      Float_t in_ptErr
     ) : pt(in_pt) ,
    eta(in_eta) ,
    phi(in_phi) ,
    charge(in_charge) ,
    ptErr(in_ptErr) {};
    virtual ~MuonFit(){};

    ClassDef(MuonFit,1)
  };

  class triMuon_cand{
  public:

	 Double_t mu_pt_max;
    Double_t mu_pt_min;
    Double_t mu_eta_max;
    Double_t mu_eta_min;

    Double_t dR_mu1mu2;
    Double_t dR_mu2mu3;
    Double_t dR_mu3mu1;

    Double_t dz_mu1mu2;
    Double_t dz_mu2mu3;
    Double_t dz_mu3mu1;

    Double_t M3mu;
    Double_t M_mu1mu2;
    Double_t M_mu2mu3;
    Double_t M_mu3mu1;

    std::array<Double_t, 3> mu_p;
    std::array<Double_t, 3> mu_pt;
    std::array<Double_t, 3> mu_pterr;
    std::array<Double_t, 3> mu_bestPt;
    std::array<Double_t, 3> mu_eta;
    std::array<Double_t, 3> mu_phi;
    std::array<Double_t, 3> mu_charge;
    std::array<Double_t, 3> vx;
    std::array<Double_t, 3> vy;
    std::array<Double_t, 3> vz;

    std::array<Double_t, 3> mu_uSta;
    std::array<Double_t, 3> mu_trkKink;
    std::array<Double_t, 3> mu_log_glbKink;
    std::array<Double_t, 3> mu_trkRelChi2;
    std::array<Double_t, 3> mu_staRelChi2;
    std::array<Double_t, 3> mu_Chi2LP;
    std::array<Double_t, 3> mu_Chi2LM;
    std::array<Double_t, 3> mu_localDist;
    std::array<Double_t, 3> mu_glbDEP;
    std::array<Double_t, 3> mu_tightMatch;
    std::array<Double_t, 3> mu_glbTrkProb;
    std::array<Double_t, 3> mu_calEM;
    std::array<Double_t, 3> mu_calEMS9;
    std::array<Double_t, 3> mu_calEMS25;
    std::array<Double_t, 3> mu_calHad;
    std::array<Double_t, 3> mu_calHadS9;
    std::array<Int_t, 3> mu_nOMS;
    std::array<Int_t, 3> mu_nOM;
    std::array<Double_t, 3> mu_comp2d;
    std::array<Double_t, 3> mu_calocomp;
    std::array<Double_t, 3> mu_segcomp;
	 
	 // Inner track information
    std::array<Double_t, 3> mu_inTrk_p;
    std::array<Double_t, 3> mu_inTrk_pt;
    std::array<Double_t, 3> mu_inTrk_eta;
    std::array<Double_t, 3> mu_inTrk_phi;
    std::array<Double_t, 3> mu_inTrk_nC2;
    std::array<Int_t, 3> mu_inTrk_trkLayWithMeas;
    std::array<Int_t, 3> mu_inTrk_pixLayWithMeas;
    std::array<Int_t, 3> mu_inTrk_nHitsTracker;
    std::array<Int_t, 3> mu_inTrk_nHitsPixel;
    std::array<Int_t, 3> mu_inTrk_nLostTrkHits;
    std::array<Int_t, 3> mu_inTrk_nLostTrkHits_in;
    std::array<Int_t, 3> mu_inTrk_nLostTrkHits_out;
    std::array<Double_t, 3> mu_trkHP;
    
	 // Muon outer track info
    std::array<Double_t, 3> mu_oTrk_p; 
    std::array<Double_t, 3> mu_oTrk_pt; 
    std::array<Double_t, 3> mu_oTrk_eta; 
    std::array<Double_t, 3> mu_oTrk_phi; 
	 std::array<Int_t, 3> mu_ioTrk_nHits;

    std::array<Double_t, 3> mu_oTrk_nC2; 
    std::array<Double_t, 3> mu_oTrk_MSWVH; 
    std::array<Double_t, 3> mu_oTrk_qprod; 

    // Muon global track info
    std::array<Double_t, 3> mu_glb_p; 
    std::array<Double_t, 3> mu_glb_eta; 
    std::array<Double_t, 3> mu_glb_phi; 
    std::array<Double_t, 3> mu_glb_nC2; 
    std::array<Double_t, 3> mu_glb_nOVMH;
    std::array<Int_t, 3>  mu_glb_nHits;
    
	 // Muon ID
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
    std::array<bool, 3> mu_hasInnerTrack;
    std::array<bool, 3> mu_hasTunePTrack;
    std::array<bool, 3> mu_hasPickyTrack;
    std::array<bool, 3> mu_hasDytTrack;
    std::array<bool, 3> mu_hasTpfmsTrack;


    // Detector based isolation
    std::array<Double_t, 3> mu_IsoR03_sumPt;
    std::array<Double_t, 3> mu_IsoR03_nTrks;
    std::array<Double_t, 3> mu_IsoR03_emEt;
    std::array<Double_t, 3> mu_IsoR03_hadEt;
    std::array<Double_t, 3> mu_IsoR03_emVetoEt;
    std::array<Double_t, 3> mu_IsoR03_hadVetoEt;
	 
	 // PF Isolation
    std::array<Double_t, 3> chargedHadronIso;
    std::array<Double_t, 3> chargedHadronIsoPU;
    std::array<Double_t, 3> photonIso;
    std::array<Double_t, 3> neutralHadronIso;
	 
	 std::array<Double_t, 3> isoPflow04; // PF isolation in dR<0.4 cone dBeta
    std::array<Double_t, 3> isoPflow03; // PF isolation in dR<0.3 cone dBeta
	 
	 std::array<Double_t, 3> dxy;   // signed transverse distance to primary vertex [cm]
    std::array<Double_t, 3> dz;  // signed longitudinal distance to primary vertex at min. transv. distance [cm]
    std::array<Double_t, 3> edxy;  // uncertainty on dxy [cm]
    std::array<Double_t, 3> edz;   // uncertainty on dz [cm]
    std::array<Double_t, 3> dxybs;  // signed transverse distance to beamspot [cm]
    std::array<Double_t, 3> dzbs;  // signed longitudinal distance to beamspot [cm]
	 std::array<Double_t, 3> mu_dxyBest; // signed tranverse distance to primary vertex [cm] using best track
	 std::array<Double_t, 3> mu_dzBest; // signed longitudnal distance to primary vertex [cm] using best track
	 std::array<Double_t, 3> dxyInner;  
    std::array<Double_t, 3> dzInner; 

    // Tight muon
    std::array<bool, 3> mu_TMOST;
    std::array<bool, 3> mu_TMOSAT;
    std::array<bool, 3> mu_TMLST;
    std::array<bool, 3> mu_TMLSAT;
    std::array<bool, 3> mu_TMLSOLPT;
    std::array<bool, 3> mu_TMSOBLPT;

    std::array<Double_t, 3> mu_tau_dR;
    std::array<Double_t, 3> mu_timAtIpInOutErr;

    // Muon time 
    std::array<Double_t, 3> muonTimeDof; 
    std::array<Double_t, 3> muonTime; 
    std::array<Double_t, 3> muonTimeErr;

    // Muon time 
    std::array<Double_t, 3> muonRpcTimeDof; 
    std::array<Double_t, 3> muonRpcTime; 
    std::array<Double_t, 3> muonRpcTimeErr;

    std::vector<tau23mu::HitInfo> hits;
    std::vector<tau23mu::ChambMatch> matches;
    std::vector<tau23mu::MuonFit> fits;

    triMuon_cand(){};
    virtual ~triMuon_cand(){};
	 
	 // reset all the variables 
    virtual void reset(){

   // for (size_t i=0; i<3; i++){
   // }
    mu_pt_max = 0;
    mu_eta_max = 0;
    mu_pt_min = 0;
    mu_eta_min = 0;
    }

    inline Float_t fitPt( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).pt;
    };

    inline Float_t fitEta( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).eta;
    };

    inline Float_t fitPhi( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).phi;
    };

    inline Int_t fitCharge( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).charge;
    };

    inline Float_t fitPtErr( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).ptErr;
    };

    ClassDef(triMuon_cand,4)
  };

  class diMuonTrk_cand{
  public:

	 Double_t mu_pt_max;
    Double_t mu_pt_min;
    Double_t mu_eta_max;
    Double_t mu_eta_min;

    Double_t dR_mu1mu2;
    Double_t dR_mu2mu3;
    Double_t dR_mu3mu1;

    Double_t dz_mu1mu2;
    Double_t dz_mu2mu3;
    Double_t dz_mu3mu1;

    Double_t M3mu;
    Double_t M_mu1mu2;
    Double_t M_mu2mu3;
    Double_t M_mu3mu1;

    std::array<Double_t, 2> mu_p;
    std::array<Double_t, 2> mu_pt;
    std::array<Double_t, 2> mu_pterr;
    std::array<Double_t, 2> mu_bestPt;
    std::array<Double_t, 2> mu_eta;
    std::array<Double_t, 2> mu_phi;
    std::array<Double_t, 2> mu_charge;
    std::array<Double_t, 2> vx;
    std::array<Double_t, 2> vy;
    std::array<Double_t, 2> vz;

    std::array<Double_t, 2> mu_uSta;
    std::array<Double_t, 2> mu_trkKink;
    std::array<Double_t, 2> mu_log_glbKink;
    std::array<Double_t, 2> mu_trkRelChi2;
    std::array<Double_t, 2> mu_staRelChi2;
    std::array<Double_t, 2> mu_Chi2LP;
    std::array<Double_t, 2> mu_Chi2LM;
    std::array<Double_t, 2> mu_localDist;
    std::array<Double_t, 2> mu_glbDEP;
    std::array<Double_t, 2> mu_tightMatch;
    std::array<Double_t, 2> mu_glbTrkProb;
    std::array<Double_t, 2> mu_calEM;
    std::array<Double_t, 2> mu_calEMS9;
    std::array<Double_t, 2> mu_calEMS25;
    std::array<Double_t, 2> mu_calHad;
    std::array<Double_t, 2> mu_calHadS9;
    std::array<Int_t, 2> mu_nOMS;
    std::array<Int_t, 2> mu_nOM;
    std::array<Double_t, 2> mu_comp2d;
    std::array<Double_t, 2> mu_calocomp;
    std::array<Double_t, 2> mu_segcomp;
	 
	 // Inner track information
    std::array<Double_t, 2> mu_inTrk_p;
    std::array<Double_t, 2> mu_inTrk_pt;
    std::array<Double_t, 2> mu_inTrk_eta;
    std::array<Double_t, 2> mu_inTrk_phi;
    std::array<Double_t, 2> mu_inTrk_nC2;
    std::array<Int_t, 2> mu_inTrk_trkLayWithMeas;
    std::array<Int_t, 2> mu_inTrk_pixLayWithMeas;
    std::array<Int_t, 2> mu_inTrk_nHitsTracker;
    std::array<Int_t, 2> mu_inTrk_nHitsPixel;
    std::array<Int_t, 2> mu_inTrk_nLostTrkHits;
    std::array<Int_t, 2> mu_inTrk_nLostTrkHits_in;
    std::array<Int_t, 2> mu_inTrk_nLostTrkHits_out;
    std::array<Double_t, 2> mu_trkHP;
    
	 // Muon outer track info
    std::array<Double_t, 2> mu_oTrk_p; 
    std::array<Double_t, 2> mu_oTrk_pt; 
    std::array<Double_t, 2> mu_oTrk_eta; 
    std::array<Double_t, 2> mu_oTrk_phi; 
	 std::array<Int_t, 2> mu_ioTrk_nHits;

    std::array<Double_t, 2> mu_oTrk_nC2; 
    std::array<Double_t, 2> mu_oTrk_MSWVH; 
    std::array<Double_t, 2> mu_oTrk_qprod; 

    // Muon global track info
    std::array<Double_t, 2> mu_glb_p; 
    std::array<Double_t, 2> mu_glb_eta; 
    std::array<Double_t, 2> mu_glb_phi; 
    std::array<Double_t, 2> mu_glb_nC2; 
    std::array<Double_t, 2> mu_glb_nOVMH;
    std::array<Int_t, 2>  mu_glb_nHits;
    
	 // Muon ID
	 std::array<bool, 2> mu_isPF;
    std::array<bool, 2> mu_isRPC;
    std::array<bool, 2> mu_isStandAlone;
    std::array<bool, 2> mu_isTrackerArb;
    std::array<bool, 2> mu_isTracker;
    std::array<bool, 2> mu_isGlobal;
    std::array<bool, 2> mu_isSoft;
    std::array<bool, 2> mu_isTight;
    std::array<bool, 2> mu_isLoose;
    std::array<bool, 2> mu_isMedium;
    std::array<bool, 2> mu_isHighPt;
    std::array<bool, 2> mu_hasInnerTrack;
    std::array<bool, 2> mu_hasTunePTrack;
    std::array<bool, 2> mu_hasPickyTrack;
    std::array<bool, 2> mu_hasDytTrack;
    std::array<bool, 2> mu_hasTpfmsTrack;


    // Detector based isolation
    std::array<Double_t, 2> mu_IsoR03_sumPt;
    std::array<Double_t, 2> mu_IsoR03_nTrks;
    std::array<Double_t, 2> mu_IsoR03_emEt;
    std::array<Double_t, 2> mu_IsoR03_hadEt;
    std::array<Double_t, 2> mu_IsoR03_emVetoEt;
    std::array<Double_t, 2> mu_IsoR03_hadVetoEt;
	 
	 // PF Isolation
    std::array<Double_t, 2> chargedHadronIso;
    std::array<Double_t, 2> chargedHadronIsoPU;
    std::array<Double_t, 2> photonIso;
    std::array<Double_t, 2> neutralHadronIso;
	 
	 std::array<Double_t, 2> isoPflow04; // PF isolation in dR<0.4 cone dBeta
    std::array<Double_t, 2> isoPflow03; // PF isolation in dR<0.3 cone dBeta
	 
	 std::array<Double_t, 2> dxy;   // signed transverse distance to primary vertex [cm]
    std::array<Double_t, 2> dz;  // signed longitudinal distance to primary vertex at min. transv. distance [cm]
    std::array<Double_t, 2> edxy;  // uncertainty on dxy [cm]
    std::array<Double_t, 2> edz;   // uncertainty on dz [cm]
    std::array<Double_t, 2> dxybs;  // signed transverse distance to beamspot [cm]
    std::array<Double_t, 2> dzbs;  // signed longitudinal distance to beamspot [cm]
	 std::array<Double_t, 2> mu_dxyBest; // signed tranverse distance to primary vertex [cm] using best track
	 std::array<Double_t, 2> mu_dzBest; // signed longitudnal distance to primary vertex [cm] using best track
	 std::array<Double_t, 2> dxyInner;  
    std::array<Double_t, 2> dzInner; 

    // Tight muon
    std::array<bool, 2> mu_TMOST;
    std::array<bool, 2> mu_TMOSAT;
    std::array<bool, 2> mu_TMLST;
    std::array<bool, 2> mu_TMLSAT;
    std::array<bool, 2> mu_TMLSOLPT;
    std::array<bool, 2> mu_TMSOBLPT;

    std::array<Double_t, 2> mu_tau_dR;
    std::array<Double_t, 2> mu_timAtIpInOutErr;

    // Muon time 
    std::array<Double_t, 2> muonTimeDof; 
    std::array<Double_t, 2> muonTime; 
    std::array<Double_t, 2> muonTimeErr;

    // Muon time 
    std::array<Double_t, 2> muonRpcTimeDof; 
    std::array<Double_t, 2> muonRpcTime; 
    std::array<Double_t, 2> muonRpcTimeErr;

    std::vector<tau23mu::HitInfo> hits;
    std::vector<tau23mu::ChambMatch> matches;
    std::vector<tau23mu::MuonFit> fits;

    diMuonTrk_cand(){};
    virtual ~diMuonTrk_cand(){};
	 
	 // reset all the variables 
    virtual void reset(){

   // for (size_t i=0; i<3; i++){
   // }
    mu_pt_max = 0;
    mu_eta_max = 0;
    mu_pt_min = 0;
    mu_eta_min = 0;
    }

    inline Float_t fitPt( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).pt;
    };

    inline Float_t fitEta( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).eta;
    };

    inline Float_t fitPhi( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).phi;
    };

    inline Int_t fitCharge( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).charge;
    };

    inline Float_t fitPtErr( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).ptErr;
    };

    ClassDef(diMuonTrk_cand,4)
  };


  class GenInfo {
  public:
    Float_t trueNumberOfInteractions; // Number of simultaneous interactions generated (before poissonian ev by ev smearing)
    Int_t actualNumberOfInteractions; // Number of simultaneous interactions generated (after poissonian ev by ev smearing)
    Float_t genWeight;
    GenInfo(){};
    virtual ~GenInfo(){};

    ClassDef(GenInfo,1)
  };

  class GenParticle {
  public:

    GenParticle(){};
    virtual ~GenParticle(){};

    Int_t pdgId;  // PDG identifier
    Int_t status; // MC status
    Float_t energy; // energy [GeV]
    Float_t pt; // pt [GeV]
    Float_t eta; // eta
    Float_t phi; // phi
    Float_t vx; // x coordinate of production vertex [cm]
    Float_t vy; // y coordinate of production vertex [cm]
    Float_t vz; // z coordinate of production vertex [cm]
    std::vector<Int_t> mothers; // vector of indices of mothers
    std::vector<Int_t> daughters; //vector of indices of daughters
    std::vector<bool>  flags; // vector of flags, in the same order of
    //  of "DataFormats/HepMCCandidate/interface/GenStatusFlag.h"

  private:

    ClassDef(GenParticle,1)
  };

  class METs {
  public:
    Float_t pfMet; // raw PF MET [GeV]
    Float_t pfChMet; // raw PF charged MET [GeV]
    Float_t caloMet; // raw Calo MET [GeV]

    METs(){};
    virtual ~METs(){};

    ClassDef(METs,1)
  };


  class Muon {
  public:

    Float_t pt;  // pt [GeV] 
    Float_t eta; // eta
    Float_t phi; // phi

    Int_t charge;  // charge

    Int_t isGlobal;
    Int_t isTracker;
    Int_t isTrackerArb;
    Int_t isRPC;
    Int_t isStandAlone;
    Int_t isPF;

    Int_t isSoft;
    Int_t isLoose;
    Int_t isTight;
    Int_t isMedium;
    Int_t isHighPt;

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

    Float_t dxy;   // signed transverse distance to primary vertex [cm]
    Float_t dz;    // signed longitudinal distance to primary vertex at min. transv. distance [cm]
    Float_t edxy;  // uncertainty on dxy [cm]
    Float_t edz;   // uncertainty on dz [cm]
    Float_t dxybs;   // signed transverse distance to beamspot [cm]
    Float_t dzbs;  // signed longitudinal distance to beamspot [cm]
    Float_t vx;
    Float_t vy;
    Float_t vz;

    Int_t nHitsGlobal;
    Int_t nHitsTracker;
    Int_t nHitsStandAlone; 

    // Variables for ID 
    //  - General (Tight, HighPt, Soft) 
    Float_t glbNormChi2; 
    Float_t trkNormChi2; 
    Int_t trkMuonMatchedStations; 
    Int_t glbMuonValidHits; 
    Int_t trkPixelValidHits; 
    Int_t trkPixelLayersWithMeas; 
    Int_t trkTrackerLayersWithMeas; 

    //  - HighPt 
    Float_t bestMuPtErr; 

    //  - Medium 
    Float_t trkValidHitFrac; 
    Float_t trkStaChi2; 
    Float_t trkKink; 
    Float_t muSegmComp; 

    //  - Soft 
    Int_t isTrkMuOST; 
    Int_t isTrkHP; 

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

    Muon(){};
    virtual ~Muon(){};

    inline Float_t fitPt( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).pt;
    };

    inline Float_t fitEta( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).eta;
    };

    inline Float_t fitPhi( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).phi;
    };

    inline Int_t fitCharge( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).charge;
    };

    inline Float_t fitPtErr( const tau23mu::MuonFitType type ) 
    {
    return fits.at(type).ptErr;
    };

    ClassDef(Muon,4)
  };

  class Track {
  public:

    Float_t pt;
    Float_t eta;
    Float_t phi;
    Float_t nChi2;
    Float_t charge;
    UInt_t nValidHits;
    Float_t dxy;
    Float_t dz;
    Float_t dxyError;
    Float_t dzError;
    Float_t vx;
    Float_t vy;
    Float_t vz;

    Track(){};
    virtual ~Track(){};

    ClassDef(Track,1)
  };

  class HLTObject {
  public:

    std::string filterTag; // name of filter passed by the object
    Float_t pt;    // pt of the object passing the filter [GeV]
    Float_t eta;     // eta of the object passing the filter
    Float_t phi;     // phi of the object passing the filter

    HLTObject(){};
    virtual ~HLTObject(){};

    ClassDef(HLTObject,1)

  };

  class L1Muon {
  public:

    Float_t pt;  // pt [GeV]
    Float_t eta; // eta
    Float_t phi; // phi
    Int_t charge; //charge (0 if invalid)

    Int_t quality;
    Int_t bx;

    Int_t tfIndex;

    L1Muon(){};
    virtual ~L1Muon(){};

    ClassDef(L1Muon,1)

  };

  class HLT {
  public:
    std::vector<std::string> triggers; // vector of strings with HLT paths
    std::vector<tau23mu::HLTObject> objects;  // vector of hlt objects assing filters

    HLT(){};
    virtual ~HLT(){};
    bool match( const std::string & path ) {
    if (  std::find (  triggers.begin(), triggers.end(), path ) != triggers.end() )
      return true;

    return false;
    }

    bool find( const std::string & path ) {
    for ( std::vector<std::string>::const_iterator it = triggers.begin(); it != triggers.end(); ++it ) {
      if ( it->find ( path ) != std::string::npos ) return true;
    }
    return false;
    }

    ClassDef(HLT,1)

  };

  class EventId {
  public:

    Int_t runNumber;      // run number
    Int_t luminosityBlockNumber;  // luminosity block number
    Int_t eventNumber;    // event number
    Int_t nMuons;       // number of good muons in the event
    std::vector<Float_t> maxPTs;  // max PT for each good muon in the event

    EventId(){};
    virtual ~EventId(){};

    ClassDef(EventId,2)
  };


  class Event {
  public:

    Int_t runNumber;     // run number
    Int_t luminosityBlockNumber; // luminosity block number
    Int_t eventNumber;     // event number

    Int_t bxId;      // bunch crossing number
    unsigned long long orbit;  // orbit number
    Float_t instLumi;    // inst lumi from scalers [10E30]

    Int_t nVtx;        // number of valid reconstructed primary vertices
    Int_t nTrk;    // number of tracks
	 Int_t n3mu; // number of 3mu candidates
	 Int_t n2muTrk; // number of 2mu + trk candidates
    Float_t primaryVertex[3];    // 3d coordinates of PV [cm]
    Float_t cov_primaryVertex[3][3]; // 3x3 covariance matrix of PV estimation [cm*cm]

    std::vector<tau23mu::GenInfo> genInfos;    // venctor of genInfos; size=0 in data
    std::vector<tau23mu::GenParticle> genParticles; // venctor of genParticles size=0 in data
    std::vector<tau23mu::Muon> muons; // vector of muons
	 std::vector<tau23mu::triMuon_cand> triMuon_cand_coll;
	 std::vector<tau23mu::diMuonTrk_cand> diMuonTrk_cand_coll;
    std::vector<tau23mu::Track> tracks; // vector of tracks
    tau23mu::METs mets;  // vector of different MET definitions 
    tau23mu::HLT hlt;       // HLT objects
    std::vector <tau23mu::L1Muon> l1muons; //vector with the L1 muon candidates



    Event(){};
    virtual ~Event(){};

    ClassDef(Event,6)
  };

}
#endif
