//// Defination of branches use dfor Tau23Mu analysis
//	  
//   Original Package:  MyAnalyzer/dsTreeMaker
//   Original Author:	Jian Wang
//   Created:  Tue, 21 Apr 2015 08:27:25 GMT

class AnalysisTreeContainer {
	public: 

		AnalysisTreeContainer (){
			edm::Service<TFileService> fs;
			tree = fs->Make<TTree>("DsTau23MuTree","DSTAU23MUTREE");
		}

	void writeTree(){
		tree->Write();
	}

	void fillTree(){
		tree->Fill();
	}
	
	TTree* tree;
};

class fillAnalysisTree(){
	public:
			
}

class muonKinematics(){
	public:
		
}

class bTagInfo(){
	public:
}

class jetInfo(){
	public:
}


