#include <iostream> 
#include <math.h>
#include "tdrstyle.C"
#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#define nmax 100

using namespace std;

void alpha_distribution(){
	
	setTDRStyle();
	
	TFile *data_file = new TFile("tree_data_2016.root","READ");
	TFile *signal_file = new TFile("tree_mc_signal.root","READ");
	TDirectory *data_dir = (TDirectory*)data_file->Get("demo");
	TDirectory *signal_dir = (TDirectory*)signal_file->Get("demo");
	TTree *data_tree = (TTree*)data_dir->Get("data");
	TTree *signal_tree = (TTree*)signal_dir->Get("data");
	
	TCanvas *c1 = new TCanvas("c1","Signal/Background",1200,1200);
	
	TH2F *h2_data = new TH2F("h2_data","Background",100,0,5,100,0,5);	
	TH2F *h2_signal = new TH2F("h2_signal","Signal",100,0,0.03,100,0,10);	
	TH1F *h_data = new TH1F("h_data","Background",100,0,0.2);
	TH1F *h_signal = new TH1F("h_signal","Signal",100,0,0.2);
	
	Int_t n_reco, n_sv;
	Double_t hlt_doublemu3_tau3mu, fv_nC;
	Double_t fv_cosdphi3D, fv_d3D, m3mu_reco, eta_reco[3], phi_reco[3], p_reco[3];	
	
	TLegend *leg = new TLegend(0.7,0.4,0.89,0.45);
	
	data_tree->SetBranchAddress("n_reco",&n_reco);
	data_tree->SetBranchAddress("hlt_doublemu3_tau3mu",&hlt_doublemu3_tau3mu);
	data_tree->SetBranchAddress("fv_cosdphi3D",&fv_cosdphi3D);
	data_tree->SetBranchAddress("fv_d3D",&fv_d3D);
	data_tree->SetBranchAddress("m3mu_reco",&m3mu_reco);
	data_tree->SetBranchAddress("eta_reco",eta_reco);
	data_tree->SetBranchAddress("fv_nC",&fv_nC);
	data_tree->SetBranchAddress("phi_reco",phi_reco);
	data_tree->SetBranchAddress("p_reco",p_reco);
	
		for (Int_t i=0; i<(data_tree->GetEntries()); ++i){
			data_tree->GetEntry(i);
			if (n_reco>=3 && hlt_doublemu3_tau3mu==1) h_data->Fill(acos(fv_cosdphi3D));//h_data->Fill(acos(fv_cosdphi3D));//h_data->Fill(fv_nC);{h2_data->Fill(eta_reco[0],phi_reco[0]);}	
	}
	h_data->Scale(1/(h_data->GetEntries()));
	h_data->Draw();
	
	signal_tree->SetBranchAddress("n_reco",&n_reco);
	signal_tree->SetBranchAddress("hlt_doublemu3_tau3mu",&hlt_doublemu3_tau3mu);
	signal_tree->SetBranchAddress("fv_cosdphi3D",&fv_cosdphi3D);
	signal_tree->SetBranchAddress("fv_d3D",&fv_d3D);
	signal_tree->SetBranchAddress("fv_nC",&fv_nC);
	signal_tree->SetBranchAddress("m3mu_reco",&m3mu_reco);
	
		for (Int_t i=0; i<(signal_tree->GetEntries()); ++i){
			signal_tree->GetEntry(i);
			if (n_reco>=3 && hlt_doublemu3_tau3mu==1) h_signal->Fill(acos(fv_cosdphi3D));//h_signal->Fill(fv_nC);
		}	
	//h2_data->Draw("COLZ");
	//h2_data->GetXaxis()->SetTitle("fv_nC");
	//h2_data->GetYaxis()->SetTitle("fv_d3D");
	//c1->SaveAs("fv_nc_vs_d3D.pdf");
	/*c1->Clear();
	h2_signal->Draw("COLZ");
	h2_signal->GetXaxis()->SetTitle("log(1+#alpha)");
	h2_signal->GetYaxis()->SetTitle("#chi^2");
	c1->SaveAs("alpha_singal.pdf");
	*/
	
	h_signal->Scale(1/(h_signal->GetEntries()));
	h_signal->SetLineColor(kRed);
	h_signal->Draw("SAMES");
	gPad->Update();
	
	h_data->SetMaximum(0.15);
	TPaveStats *st_data = (TPaveStats*)h_data->FindObject("stats");
	st_data->SetX1NDC(0.7);
	st_data->SetX2NDC(0.9);
	st_data->SetY1NDC(0.7);
	st_data->SetY2NDC(0.9);
	
	TPaveStats *st_signal = (TPaveStats*)h_signal->FindObject("stats");
	st_signal->SetX1NDC(0.7);
	st_signal->SetX2NDC(0.9);
	st_signal->SetY1NDC(0.5);
	st_signal->SetY2NDC(0.7);
	
	leg->AddEntry("h_data");
	leg->AddEntry("h_signal");
	leg->Draw("same");
	c1->SaveAs("alpha.pdf");
	//c1->SaveAs("p_reco.pdf");
	//c1->SaveAs("Normalized_ChiSquared.pdf");
	//c1->SaveAs("vertex_displacement.pdf");
	
	}
