#include <iostream>

using namespace edm;
using namespace std;

pass ESHandle??

TransientVertex fitVertex(TrackRef _trk1, TrackRef _trk2, ,TrackRef _trk3){

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
      
      if(!fv.isValid()) { return; cout<<"Vertex Fit unvalid!"<<endl; 
		}

		return fv;
}
