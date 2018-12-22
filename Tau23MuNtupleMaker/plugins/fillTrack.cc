void MuonPogTreeProducer::fillTracks(const edm::Handle<std::vector<reco::Track> >& tracks,
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
