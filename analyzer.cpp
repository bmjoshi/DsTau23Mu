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

void analyzer(){

	setTDRStyle();
	
	TFile *data_file = new TFile("tree_data.root","READ");
	TFile *signal_file = new TFile("tree_signal.root","READ");
	TDirectory *data_dir = (TDirectory*)data_file->Get("demo");
	TDirectory *signal_dir = (TDirectory*)signal_file->Get("demo");
	TTree *data_tree = (TTree*)data_dir->Get("data");
	TTree *signal_tree = (TTree*)signal_dir->Get("data");
	
	TCanvas *c1 = new TCanvas("c1","Signal/ Background",1200,1200);
	
	TH1F *h_data = new TH1F("h_data","Background",100,1.6,2.0);
	TH1F *h_signal = new TH1F("h_signal","Signal",100,1.6,2.0);
	TH1F *h_m3mu = new TH1F("h_m3mu", "",100,1.6,2.);
	
	TLegend *leg = new TLegend(0.66,0.66,0.81,0.74);
	TLegend *leg_m3mu = new TLegend(0.66,0.76,0.81,0.84);
	Int_t n_reco;
	Double_t nOMS_min, p_reco[3], pt_reco[3], sv_cosdphi3D, hlt_doublemu3_tau3mu, m3mu_reco, ntrk_tau05, dr12_reco, dr23_reco, dr31_reco, lumiN, eta_reco[3], phi_reco[3], pTx, pTy, pz, m_inv, mT_inv, e, eT, Iso_nTr;
	Double_t mu_mass = 0.10586;
	Float_t IntLumi;

	TObjArray* arr = (TObjArray*)data_tree->GetListOfBranches();
	for (int i=0; i<arr->GetSize(); i++){
		TString str = arr->At(i)->GetName();
		Double_t max = data_tree->GetMaximum(str);
		Double_t min = data_tree->GetMinimum(str);
		Int_t entries = data_tree->GetEntries(str);
		cout<<str<<" "<<min<<" "<<max<<endl;
		} 
	//data_tree->Draw("fv_tC>>h_data","n_reco>=3 && hlt_doublemu3_tau3mu==1");
	//signal_tree->Draw("fv_tC>>h_signal","n_reco>=3 && hlt_doublemu3_tau3mu==1","SAME");
	
	data_tree->SetBranchAddress("lumiN",&lumiN);
	data_tree->SetBranchAddress("m3mu_reco",&m3mu_reco);
	data_tree->SetBranchAddress("n_reco",&n_reco);
	data_tree->SetBranchAddress("p_reco",p_reco);
	data_tree->SetBranchAddress("pt_reco",pt_reco);
	data_tree->SetBranchAddress("eta_reco",eta_reco);
	data_tree->SetBranchAddress("phi_reco",phi_reco);
	data_tree->SetBranchAddress("hlt_doublemu3_tau3mu",&hlt_doublemu3_tau3mu);
	data_tree->SetBranchAddress("dr12_reco",&dr12_reco);
	data_tree->SetBranchAddress("dr23_reco",&dr23_reco);
	data_tree->SetBranchAddress("dr31_reco",&dr31_reco);
	data_tree->SetBranchAddress("Iso_nTr",&Iso_nTr);
	data_tree->SetBranchAddress("nOMS_min",&nOMS_min);
	data_tree->SetBranchAddress("sv_cosdphi3D",sv_cosdphi3D);
/*
	// Calculate Integrated luminosity
	// Plot signal against background
	for (Int_t i=0; i<(data_tree->GetEntries()); i++) {
		data_tree->GetEntry(i);
		if (m3mu_reco>=1.6 && m3mu_reco<=2.0 && n_reco>=3 && hlt_doublemu3_tau3mu==1) h_data->Fill(m3mu_reco);
		}
	signal_tree->SetBranchAddress("m3mu_reco",&m3mu_reco);
	signal_tree->SetBranchAddress("n_reco",&n_reco);
	signal_tree->SetBranchAddress("hlt_doublemu3_tau3mu",&hlt_doublemu3_tau3mu);
	for (Int_t i=0; i<(signal_tree->GetEntries()); i++) {
		signal_tree->GetEntry(i);
		if (m3mu_reco>=1.6 && m3mu_reco<=2.0 && n_reco>=3 && hlt_doublemu3_tau3mu==1) h_signal->Fill(m3mu_reco);
		}
*/
/*	for (Int_t i=0; i<(data_tree->GetEntries()); i++) {
		data_tree->GetEntry(i);
		if (n_reco>=3 && hlt_doublemu3_tau3mu==1){
			//h_data->Fill(sv_cosdphi3D);
			}
		}
	h_data->Scale(1/(h_data->GetEntries()));

	signal_tree->SetBranchAddress("lumiN",&lumiN);
	signal_tree->SetBranchAddress("m3mu_reco",&m3mu_reco);
	signal_tree->SetBranchAddress("n_reco",&n_reco);
	signal_tree->SetBranchAddress("p_reco",p_reco);
	signal_tree->SetBranchAddress("pt_reco",pt_reco);
	signal_tree->SetBranchAddress("eta_reco",eta_reco);
	signal_tree->SetBranchAddress("phi_reco",phi_reco);
	signal_tree->SetBranchAddress("hlt_doublemu3_tau3mu",&hlt_doublemu3_tau3mu);
	signal_tree->SetBranchAddress("dr12_reco",&dr12_reco);
	signal_tree->SetBranchAddress("dr23_reco",&dr23_reco);
	signal_tree->SetBranchAddress("dr31_reco",&dr31_reco);
	signal_tree->SetBranchAddress("Iso_nTr",&Iso_nTr);
	signal_tree->SetBranchAddress("ntrk_tau05",&ntrk_tau05);
	signal_tree->SetBranchAddress("nOMS_min",&nOMS_min);
	signal_tree->SetBranchAddress("sv_cosdphi3D",sv_cosdphi3D);
	
	//std::string parLabel[] = {"lumiN", "m3mu_reco", "p_reco", "eta_reco", };
	//Double_t histMin[] = {};
	//Double_t histMax[] = {};
	
	for (Int_t i=0; i<(signal_tree->GetEntries()); i++) {
		signal_tree->GetEntry(i);
		if (n_reco>=3 && hlt_doublemu3_tau3mu==1){
			//h_signal->Fill(sv_cosdphi3D); 			
			}
		}

	h_signal->Scale(1/(h_signal->GetEntries()));	
	h_signal->Draw();
	h_signal->SetMaximum(0.42);
	h_signal->SetLineColor(kRed);
	h_signal->GetYaxis()->SetTitle("a.u.");
	h_signal->GetYaxis()->SetTitleOffset(2.0);
	h_signal->GetXaxis()->SetTitle("3#mu vertex, #chi^2 dof");
	
	//h_data->SetTitle("");
	h_data->Draw("SAME");
	h_data->SetMarkerSize(2);
	h_data->SetMarkerStyle(20);

	leg->AddEntry(h_data, "data (2016)");
	leg->AddEntry(h_signal, "MC data");
	leg->Draw("same");

	c1->Modified();
	c1->Update();
	c1->SaveAs("fv_tC_svg.png");

	//Plot other parameters
	//Plot m3mu_reco & m_inv
	
	for (Int_t i=0; i<(data_tree->GetEntries()); i++) {
		data_tree->GetEntry(i);
		if (n_reco>=3 && hlt_doublemu3_tau3mu==1){
			pz=pTx=pTy=eT=e=0;
			for (int j=0; j<3; j++){
			pz += pt_reco[j]*sinh(eta_reco[j]);
			pTx += pt_reco[j]*cos(phi_reco[j]);
			pTy += pt_reco[j]*sin(phi_reco[j]);
			eT += sqrt(pow(mu_mass,2)+pow(pt_reco[j],2));
			e += sqrt(pow(mu_mass,2)+pow(pt_reco[j]*cosh(eta_reco[j]),2));
			}
			mT_inv = sqrt(pow(eT,2)-pow(pTy,2)-pow(pTx,2));
			m_inv = sqrt(pow(e,2)-pow(pTx,2)-pow(pTy,2)-pow(pz,2));
			h_m3mu->Fill(m_inv); 
			h_data->Fill(m3mu_reco);
			}
		}
	
	h_data->SetTitle("Invariant mass#{3#mu}");
	h_data->Draw();
	h_data->SetMarkerSize(2);
	h_data->SetMarkerStyle(20);
	h_m3mu->Draw("SAME");
	h_m3mu->SetLineColor(kRed);
	h_data->GetYaxis()->SetTitle("a.u.");
	h_data->GetYaxis()->SetTitleOffset(2.0);
	h_data->GetXaxis()->SetTitle("m3#mu (GeV)");
	h_data->Scale(1/(h_data->GetEntries()));
	h_m3mu->Scale(1/(h_m3mu->GetEntries()));
	h_data->SetMaximum(0.1);
	signal_tree->Draw("m3mu_reco>>h_signal","n_reco>=3 && hlt_doublemu3_tau3mu==1","SAME");
	h_signal->SetLineColor(kGreen+3);
	h_signal->Scale(1/(h_signal->GetEntries()));
	leg_m3mu->AddEntry(h_data, "m3mu_reco");
	leg_m3mu->AddEntry(h_m3mu, "m_inv");
	leg_m3mu->AddEntry(h_signal, "m3mu_signal");
	leg_m3mu->Draw("same");

	c1->Modified();
	c1->Update();
	c1->SaveAs("invariant_mass_svb.png");
*/
}
