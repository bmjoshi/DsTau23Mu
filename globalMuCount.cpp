#include <iostream>
#include <math.h>
#include "tdrstyle.C"
#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

using namespace std;

void globalMuCount(){

	setTDRStyle();
	
	TFile *signal_file = new TFile("tree_mc_signal.root","READ");
	TDirectory *signal_dir = (TDirectory*)signal_file->Get("demo");
	TTree *signal_tree = (TTree*)signal_dir->Get("data");
	
	TCanvas *c1 = new TCanvas("c1","Signal/ Background",1200,1200);

	Int_t n_reco;
	Double_t nOMS_min, p_reco[3], pt_reco[3], ggm[3], tgm[3], sv_cosdphi3D, hlt_doublemu3_tau3mu, m3mu_reco, ntrk_tau05, dr12_reco, dr23_reco, dr31_reco, lumiN, eta_reco[3], phi_reco[3], pTx, pTy, pz, m_inv, mT_inv, e, eT, Iso_nTr;
	Double_t mu_mass = 0.10586;
	Float_t IntLumi;
	int counter=0, hlt_dm3t3m_count=0;
	int gmuCount[4];
	Int_t nentries  = signal_tree->GetEntries();
	
	signal_tree->SetBranchAddress("lumiN",&lumiN);
	signal_tree->SetBranchAddress("m3mu_reco",&m3mu_reco);
	signal_tree->SetBranchAddress("n_reco",&n_reco);
	signal_tree->SetBranchAddress("ggm_reco",ggm);
	signal_tree->SetBranchAddress("tgm_reco",tgm);
	signal_tree->SetBranchAddress("p_reco",p_reco);
	signal_tree->SetBranchAddress("pt_reco",pt_reco);
	signal_tree->SetBranchAddress("eta_reco",eta_reco);
	signal_tree->SetBranchAddress("phi_reco",phi_reco);
	signal_tree->SetBranchAddress("hlt_doublemu3_tau3mu",&hlt_doublemu3_tau3mu);
	signal_tree->SetBranchAddress("dr12_reco",&dr12_reco);
	signal_tree->SetBranchAddress("dr23_reco",&dr23_reco);
	signal_tree->SetBranchAddress("dr31_reco",&dr31_reco);
	signal_tree->SetBranchAddress("Iso_nTr",&Iso_nTr);
	signal_tree->SetBranchAddress("nOMS_min",&nOMS_min);
	signal_tree->SetBranchAddress("sv_cosdphi3D",sv_cosdphi3D);
	
	
	//for (int k=3; k<8; k++){
	hlt_dm3t3m_count=0;
	for (int j=0; j<4; j++) gmuCount[j]=0;
	for (UInt_t i=0; i<nentries; i++){
		signal_tree->GetEntry(i);
		counter=0;
		if (n_reco>=3 && hlt_doublemu3_tau3mu==1){
			hlt_dm3t3m_count++;
			for (int j=0; j<3; j++) if (tgm[j]==1) counter++;
			gmuCount[counter]++;
			}
		}
	cout<<"Number of entries: "<<nentries<<endl;
	cout<<"Number of events passing the hlt_doublemu3_tau3mu cut and n_reco>=3: "<<hlt_dm3t3m_count<<endl;
	for (int j=0; j<4; j++) cout<<gmuCount[j]<<"\t";
	cout<<endl;
	//}
}
