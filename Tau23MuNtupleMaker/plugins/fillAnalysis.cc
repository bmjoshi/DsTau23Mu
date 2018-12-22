#include <iostream>
#include "AnalysisTree.h"
#include "Utils.h"
#include "3mu_cand.h"

//// Loop over all the muons

using namespace reco;
using namespace edm;
using namespace std;

#define dz_mumu_cut 0.5
#define dR_mumu_cut 0.8

double n_3mu
vector<Int_t> goodMuIndex;
vector<Int_t> muMatchedTrk;

vector<array<Double_t, 3> > mu_pt;
vector<array<Double_t, 3> > mu_eta;
vector<array<Double_t, 3> > mu_phi;

vector<Double_t> dR_mumu_12;
vector<Double_t> dR_mumu_23;
vector<Double_t> dR_mumu_31;

vector<Double_t> dz_mumu_12;
vector<Double_t> dz_mumu_23;
vector<Double_t> dz_mumu_31;

vector<Double_t> mu_pt_max;
vector<Double_t> mu_eta_min;
vector<Double_t> mu_pt_max;
vector<Double_t> mu_eta_min;

3mu_cand tmp_3mu_cand;
2muTrk_cand tmp_2muTrk_cand;

(const edm::Handle<edm::View<reco::Muon> > & muons,
 const edm::Handle<std::vector<reco::Vertex> > & vertexes,
 const edm::Handle<reco::BeamSpot> & beamSpot)
{

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

   for (size_t i=0; i<goodMuIndex.size()-1; ++i)s{
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
          
			 auto fv = tau23mu::fitVertex(trk1, trk2, trk3);
			 double fvnC2_tmp = fv.totalChiSquared()/fv.degreesOfFreedom();

          if (fvnC2_tmp > 3muFitVtx_nC2_cut ) continue; // Store all the good candidates

          // Sort the muons by P

          vector<reco::Muons> sortedMu;

          sortedMu.push_back(maxPMu(maxPMu(mu1,mu2),mu3));
          sortedMu.push_back(minPMu(minPMu(maxPMu(mu1,mu2),maxPMu(mu2,mu3)),maxPMu(mu3,mu1)));
          sortedMu.push_back(minPMu(minPMu(mu1,mu2),mu3));

          // Fill tree with trimuon information

          tmp_mu_pt[0] = sotredMu[0].pt();
          tmp_mu_pt[1] = sotredMu[1].pt();
          tmp_mu_pt[2] = sotredMu[2].pt();

          tmp_mu_eta[0] = sotredMu[0].eta();
          tmp_mu_eta[1] = sotredMu[1].eta();
          tmp_mu_eta[2] = sotredMu[2].eta();

          tmp_mu_phi[0] = sotredMu[0].phi();
          tmp_mu_phi[1] = sotredMu[1].phi();
          tmp_mu_phi[2] = sotredMu[2].phi();

          mu_pt_max =  TMath::Max(TMath::Max(mu1.pt(),mu2.pt()),mu3.pt()); //return highest pt muon 
          mu_eta_max = TMath::Max(TMath::Max(mu1.eta(),mu2.eta()),mu3.eta()); //return highest eta muon 
          
			 mu_pt_max =  TMath::Min(TMath::Min(mu1.pt(),mu2.pt()),mu3.pt()); //return lowest pt muon 
          mu_eta_max = TMath::Min(TMath::Min(mu1.eta(),mu2.eta()),mu3.eta()); //return lowest eta muon 
          
			 // Re-assign and store dR and dz with sorted muons

          dR_mu1mu2 = deltaR(sortedMu[0].eta(), sortedMu[1].eta(), sortedMu[0].phi(), sortedMu[1].phi());
          dR_mu2mu3 = deltaR(sortedMu[1].eta(), sortedMu[2].eta(), sortedMu[1].phi(), sortedMu[2].phi());
          dR_mu1mu3 = deltaR(sortedMu[0].eta(), sortedMu[2].eta(), sortedMu[0].phi(), sortedMu[2].phi());

          dz_mu1mu2 = abs(sortedMu[0].InnerTracker()->dz(beamSpotHandle->position())-sortedMu[1].InnerTracker()->dz(beamSpotHandle->Position()));
          dz_mu2mu3 = abs(sortedMu[1].InnerTracker()->dz(beamSpotHandle->position())-sortedMu[2].InnerTracker()->dz(beamSpotHandle->Position()));
          dz_mu1mu3 = abs(sortedMu[0].InnerTracker()->dz(beamSpotHandle->position())-sortedMu[2].InnerTracker()->dz(beamSpotHandle->Position()));

          mumu_dR_12 = dR_mu1mu2;
          mumu_dR_23 = dR_mu2mu3;
          mumu_dR_31 = dR_mu1mu3;

          mumu_dz_12 = dz_mu1mu2;
          mumu_dz_23 = dz_mu2mu3;
          mumu_dz_31 = dz_mu1mu3;

			 // Compute 3 mu invariant mass
			 TLorentzVector vec_mu1, vec_mu2, vec_mu3;

			 vec_mu1.SetPtEtaPhiM(sortedMu[0].pt(), sortedMu[0].eta(), sortedMu[0].phi(), MU_MASS);
			 vec_mu2.SetPtEtaPhiM(sortedMu[1].pt(), sortedMu[1].eta(), sortedMu[1].phi(), MU_MASS);
			 vec_mu3.SetPtEtaPhiM(sortedMu[2].pt(), sortedMu[2].eta(), sortedMu[2].phi(), MU_MASS);

			 3mu_InvMass.push_back((vec_mu1+vec_mu2+vec_mu3).M());

          3muFitVtx_nC2.push_back(fvnC2_tmp);
          n_3mu++;
        }
      }

      // Find 2mu+1track
      if (do2mu_ && mu1.pt()>2.5 && mu2.pt()>2.5){
			for (; trkIt!=trkEnd; ++trkIt){
			
			const reco::Track& trk = (*trkIt);
			
			// Build tracks  and re-fit the vertex
         TrackRef trk1 = mu1.innerTrack();
         TrackRef trk2 = mu2.innerTrack();
         TrackRef trk3 = trk;

			auto fv = 

         if(!3muFitVtx.isValid()) continue;

         double fvnC2_tmp = fv.totalChiSquared()/fv.degreesOfFreedom();

         if (fvnC2_tmp > 2muTrkFitVtx_nC2_cut ) continue; // Store all the good candidates

			// Compute 2mu + trk invaiant mass
			TLorentzVector vec_mu1, vec_mu2, vec_trk;

			vec_mu1.SetPtEtaPhiM(,MU_MASS);
			vec_mu1.SetPtEtaPhiM(,MU_MASS);
			vec_trk.SetPtEtaPhiM(,PI_MASS);

			2muTrk_InvMass.push_back((vec_mu1+vec_mu2+vec_mu3).M());

			2muTrkVtx_nC2.push_back(fvnC2_tmp);
			}
      }

    }



   }
}


