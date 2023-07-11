#define G2EEAnalysis_cxx
#include "G2EEAnalysis.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cmath>

//declare histograms

TH1D* h_cutflow_ee = new TH1D("h_cutflow_ee","Cutflow",11,0.0,11);
TH1D* h_invm = new TH1D("h_invm","Invariant mass",70,0,70);
TH1D* h_pt_p = new TH1D("h_pt_p","Transverse momentum",100,0,30);
TH1D* h_pt_n = new TH1D("h_pt_n","Transverse momentum",100,0,30);
TH1D* h_eta_p = new TH1D("h_eta_p","Eta",100,-4,4);
TH1D* h_eta_n = new TH1D("h_eta_n","Eta",100,-4,4);
TH1D* h_phi_p = new TH1D("h_phi_p","Phi",100,-5,5);
TH1D* h_phi_n = new TH1D("h_phi_n","Phi",100,-5,5);
TH2D* hh_pt1_pt2 = new TH2D("hh_pt1_pt2","Transverse momentum",100,0,30,100,0,30);
TH1D* h_aco = new TH1D("h_aco","Acoplanarity",100,0,0.01);
TH1D* h_Ptll_ee = new TH1D("Pt ee","Pt ee",100,0,5);
TH2D* hh_zdc_ee = new TH2D("ZDC","ZDC",100,0,10,100,0,10);
TH1D* yy_ee = new TH1D("yy_ovf","yy",100,-5,5);
TH1D* h_f = new TH1D("h_f","F",2,0,2);
TH1D* h_f0n0n = new TH1D("h_f0n0n","F",2,0,2);
TH1D* h_fxn0n = new TH1D("h_fxn0n","F",2,0,2);
TH1D* h_fxnxn = new TH1D("h_fxnxn","F",2,0,2);

double lims[5] = {5,10,20,40,80};

TH1D* h_m_ee_all08 = new TH1D("h_m_ee_all08","",4,lims);
TH1D* h_m_ee_0n0n08 = new TH1D("h_m_ee_0n0n08","",4,lims);
TH1D* h_m_ee_xn0n08 = new TH1D("h_m_ee_xn0n08","",4,lims);
TH1D* h_m_ee_xnxn08 = new TH1D("h_m_ee_xnxn08","",4,lims);
TH1D* h_m_ee_all16 = new TH1D("h_m_ee_all16","",4,lims);
TH1D* h_m_ee_0n0n16 = new TH1D("h_m_ee_0n0n16","",4,lims);
TH1D* h_m_ee_xn0n16 = new TH1D("h_m_ee_xn0n16","",4,lims);
TH1D* h_m_ee_xnxn16 = new TH1D("h_m_ee_xnxn16","",4,lims);
TH1D* h_m_ee_all24 = new TH1D("h_m_ee_all24","",4,lims);
TH1D* h_m_ee_0n0n24 = new TH1D("h_m_ee_0n0n24","",4,lims);
TH1D* h_m_ee_xn0n24 = new TH1D("h_m_ee_xn0n24","",4,lims);
TH1D* h_m_ee_xnxn24 = new TH1D("h_m_ee_xnxn24","",4,lims);

//event selection function

Bool_t G2EEAnalysis::SelectEETrack(){

    h_cutflow_ee->Fill(0); 
    
    if(not passed_GRL) return false; //passed GRL
    h_cutflow_ee->Fill(1);
    //HLT_hi_upc_FgapAC3_mb_sptrk_exclusiveloose2_L12TAU1_VTE50
    //if(not (passed_HLT_hi_upc_FgapAC3_hi_gg_upc_L12TAU1_VTE50 || passed_HLT_hi_upc_FgapAC3_hi_gg_upc_L1TAU1_TE4_VTE200)) return false;
    if(not (passed_HLT_hi_upc_FgapAC3_mb_sptrk_exclusiveloose2_L12TAU1_VTE50)) return false; //passed trigger
    //if(not (passed_HLT_hi_upc_FgapAC3_hi_gg_upc_L12TAU1_VTE50)) return false;
    h_cutflow_ee->Fill(2);
    
    
    
    Int_t nEl0=0;
    vector<int> index;
    vector<double> e_phi;
    for(Int_t i=0; i<nElec; ++i){
        if(electron_is_LHLoose->at(i)){ 
            nEl0++;
            index.push_back(i);
            e_phi.push_back(electron_phi->at(i));
        }
    }
    
    if (nEl0!=2||index.size()!=2) return false; //2 electrons very loose
    h_cutflow_ee->Fill(3);
    
        
    
    
    Int_t nEl2=0;
    for (Int_t i=0; i<nElec; ++i){
        if((electron_eta->at(i)>=1.37 and electron_eta->at(i)<=1.52) or electron_eta->at(i)>=2.47) nEl2++;
    }
    
    if(nEl2>0) return false; //eta <1.37 or 1.52<eta<2.47
    h_cutflow_ee->Fill(4);
    
    if(mll_ee<=5.0) return false; //m>5.0
    h_cutflow_ee->Fill(5);
    
    if(nMuon!=0) return false;
    h_cutflow_ee->Fill(6); //no muons
    
 	if(track_n_LoosePrimary!=2) return false; //number of tracks loose primary: 2
    h_cutflow_ee->Fill(7);
    
    
    //pt
    Int_t nEl=0;
    for(Int_t i=0; i<nElec; ++i){
        if(electron_pt->at(i)<=2.5) nEl++;
    }
    
    if(nEl>0) return false; //pt > 3.0
    h_cutflow_ee->Fill(8);
    
    
    double delta_Phi = e_phi[0]-e_phi[1];
    if (TMath::Abs(delta_Phi) > TMath::Pi()) delta_Phi =  2*TMath::Pi() -TMath::Abs(delta_Phi); 
    double acoll = 1 - TMath::Abs(delta_Phi)/TMath::Pi();

    
    if(electron_charge->at(index[0])*electron_charge->at(index[1])>0) return false;
    h_cutflow_ee->Fill(9); //opposite charges
    //przeciwne znaki
 
 
 	if(Ptll_ee>=2.0) return false; //pair pt <2.0 GeV
 	h_cutflow_ee->Fill(10);
 	
  
   	//if (not (zdc_ene_a < 1e3 && zdc_ene_c < 1e3)) return false;
    
    
    //fill pt histograms
    
    double pt_p, pt_n;
    for (Int_t i=0;i<nElec;++i){
        if(electron_charge->at(i)>0){
            h_pt_p->Fill(electron_pt->at(i));
            pt_p=electron_pt->at(i);
        }    
        if(electron_charge->at(i)<0){
            h_pt_n->Fill(electron_pt->at(i));
            pt_n=electron_pt->at(i);
        }
    }
    
	
    hh_pt1_pt2->Fill(pt_p,pt_n);
    for (Int_t i=0;i<nElec;++i){
        if(electron_charge->at(i)>0) h_eta_p->Fill(electron_eta->at(i));
        if(electron_charge->at(i)<0) h_eta_n->Fill(electron_eta->at(i));
    }
    
    for (Int_t i=0;i<nElec;++i){ 
        if(electron_charge->at(i)>0) h_phi_p->Fill(electron_phi->at(i));
        if(electron_charge->at(i)<0) h_phi_n->Fill(electron_phi->at(i));
    }
    
    h_Ptll_ee->Fill(Ptll_ee);
    hh_zdc_ee->Fill(zdc_ene_a/(2.51e3),zdc_ene_c/(2.51e3));

    h_aco->Fill(acoll);
    h_invm->Fill(mll_ee);
    	
	//fill mass histograms according to rapidity
	
	h_f->Fill(0.5);
    if(electron_pt->at(0)>4 && mll_ee>10 && electron_pt->at(1)>4) h_f->Fill(1.5);
	
	if (not((zdc_ene_a < 1e3 && zdc_ene_c < 1e3))){
    
    		if(zdc_ene_a < 1e3 || zdc_ene_c < 1e3){
    			h_fxn0n->Fill(0.5);
    	 		if(electron_pt->at(0)>4 && mll_ee>10 && electron_pt->at(1)>4) h_fxn0n->Fill(1.5);
    	 	}
    		else{
    			h_fxnxn->Fill(0.5);
    	 		if(electron_pt->at(0)>4 && mll_ee>10 && electron_pt->at(1)>4) h_fxnxn->Fill(1.5);
    		} 
    	
    }else{
    	 h_f0n0n->Fill(0.5);
    	 if(electron_pt->at(0)>4 && mll_ee>10 && electron_pt->at(1)>4) h_f0n0n->Fill(1.5);
	}
	
   	if(abs(yll_ee)<0.8){
    
    	h_m_ee_all08->Fill(mll_ee);
    
    	if (not((zdc_ene_a < 1e3 && zdc_ene_c < 1e3))){
    
    		if(zdc_ene_a < 1e3 || zdc_ene_c < 1e3)
    			h_m_ee_xn0n08->Fill(mll_ee);
    		else h_m_ee_xnxn08->Fill(mll_ee);
    	
    	}else h_m_ee_0n0n08->Fill(mll_ee);
    	
    	
    	
    }else if(abs(yll_ee)<1.6){
    
    	h_m_ee_all16->Fill(mll_ee);
    
    	if (not((zdc_ene_a < 1e3 && zdc_ene_c < 1e3))){
    
    		if(zdc_ene_a < 1e3 || zdc_ene_c < 1e3)
    			h_m_ee_xn0n16->Fill(mll_ee);
    		else h_m_ee_xnxn16->Fill(mll_ee);
    	
    		
    	}else h_m_ee_0n0n16->Fill(mll_ee);
    	
    }else if(abs(yll_ee)<2.4){
     	h_m_ee_all24->Fill(mll_ee);
    
    	if (not((zdc_ene_a < 1e3 && zdc_ene_c < 1e3))){
    
    		if(zdc_ene_a < 1e3 || zdc_ene_c < 1e3)
    			h_m_ee_xn0n24->Fill(mll_ee);
    		else h_m_ee_xnxn24->Fill(mll_ee);
    		
    	
    		
    	}else h_m_ee_0n0n24->Fill(mll_ee);
    }

    yy_ee->Fill(yll_ee);
  	
  	
    return true;

}


void G2EEAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L G2TauEEAnalysis.C
//      root> G2TauEEAnalysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   //cutflow bin labels
   h_cutflow_ee->GetXaxis()->SetBinLabel(1,"All");
   h_cutflow_ee->GetXaxis()->SetBinLabel(2,"GRL");
   h_cutflow_ee->GetXaxis()->SetBinLabel(3,"Trigger");
   h_cutflow_ee->GetXaxis()->SetBinLabel(4,"2 electrons");
   h_cutflow_ee->GetXaxis()->SetBinLabel(5,"|#eta| < 2.47 without crack region");
   h_cutflow_ee->GetXaxis()->SetBinLabel(6,"inv_mass > 5 GeV");
   h_cutflow_ee->GetXaxis()->SetBinLabel(7,"No muons");
   h_cutflow_ee->GetXaxis()->SetBinLabel(8,"2 loose primary tracks");
   h_cutflow_ee->GetXaxis()->SetBinLabel(9,"Electron P_{t} > 2.5 GeV");
   h_cutflow_ee->GetXaxis()->SetBinLabel(10,"opposite charges");
   h_cutflow_ee->GetXaxis()->SetBinLabel(11,"Pair P_{t} < 2.0 GeV");
   
   
   //TH1D* h_eta_ee = new TH1D("h_eta_ee",";electron #eta; Events",-4.0,4.0,50);

   //h_eta_ee->Sumw2();


   Long64_t nbytes = 0, nb = 0;
   int passed=0;
   int n0n0=0,xn0n=0, xnxn=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (SelectEETrack()){
        passed++;      
	    if(passed%1000==0)std::cout << passed <<  ". Event passed!" << std::endl;
	    if(zdc_ene_a<=1e3 && zdc_ene_c<=1e3) n0n0++;
    	else if(zdc_ene_a>1e3 && zdc_ene_c>1e3)xnxn++;
    	else xn0n++;
    
	    for(Int_t i=0; i<nElec; i++){
	        if((electron_eta->at(i)>=1.37 and electron_eta->at(i)<=1.52) or electron_eta->at(i)>=2.47){
	            //h_eta_ee->Fill(electron_eta->at(i));
	        } 
	    }  
	 }
   }
   //print calculations
   cout << "Passed events: " << passed << endl;
   cout << "ZDC: n0n0: " << n0n0 << " xn0n: " << xn0n << " xnxn: " << xnxn << endl;
   
   //save histograms to file

   //TFile *outHistFile =TFile::Open("outputdata.root","RECREATE");
   //TFile *outHistFile =TFile::Open("output8M.root","RECREATE");
   //TFile *outHistFile =TFile::Open("outputP6M8.root","RECREATE");
   //TFile *outHistFile =TFile::Open("outputdata0n.root","RECREATE");
   //TFile *outHistFile =TFile::Open("output15M.root","RECREATE");
   //TFile *outHistFile =TFile::Open("outputP6M15.root","RECREATE");
   //TFile *outHistFile =TFile::Open("output4p5M7.root","RECREATE");
   //TFile *outHistFile =TFile::Open("output7M15.root","RECREATE");
   TFile *outHistFile =TFile::Open("output15Mv1.root","RECREATE");
   outHistFile->cd();
   h_cutflow_ee->Write();
   h_invm->Write();
   h_pt_p->Write();
   h_pt_n->Write();
   h_eta_p->Write();
   h_eta_n->Write();
   h_phi_p->Write();
   h_phi_n->Write();
   h_aco->Write();
   h_Ptll_ee->Write();
   hh_pt1_pt2->Write();
   hh_zdc_ee->Write();
   h_m_ee_0n0n08->Write();
   h_m_ee_xn0n08->Write();
   h_m_ee_xnxn08->Write();
   h_m_ee_all08->Write();
   h_m_ee_0n0n16->Write();
   h_m_ee_xn0n16->Write();
   h_m_ee_xnxn16->Write();
   h_m_ee_all16->Write();
   h_m_ee_0n0n24->Write();
   h_m_ee_xn0n24->Write();
   h_m_ee_xnxn24->Write();
   h_m_ee_all24->Write();
   h_f0n0n->Write();
   h_fxn0n->Write();
   h_fxnxn->Write();
   h_f->Write();
   yy_ee->Write();
   outHistFile->Close();
   
}
