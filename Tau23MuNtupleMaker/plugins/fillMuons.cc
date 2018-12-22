#include <iostream>

using namespace std;
using namespace edm;

Int_t MuonPogTreeProducer::fillMuons(const edm::Handle<edm::View<reco::Muon> > & muons,
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
					 ntupleMu.trkNormChi2	        = hasInnerTrack ? mu.innerTrack()->normalizedChi2()  : -999; 
					 ntupleMu.trkMuonMatchedStations   = isTracker     ? mu.numberOfMatchedStations()       : -999; 
					 ntupleMu.glbMuonValidHits	        = isGlobal      ? mu.globalTrack()->hitPattern().numberOfValidMuonHits()       : -999; 
					 ntupleMu.trkPixelValidHits	= hasInnerTrack ? mu.innerTrack()->hitPattern().numberOfValidPixelHits()       : -999; 
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
