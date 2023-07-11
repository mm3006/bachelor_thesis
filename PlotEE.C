
void PlotEE(){

    //declare histograms
    TH1D* h_invm_all = new TH1D("h_invm_all","Invariant mass",70,0,70);
    TH1D* h_ee_pt_p_all = new TH1D("h_ee_pt_p_all","pt_p",100,0,30);
    TH1D* h_ee_pt_n_all = new TH1D("h_ee_pt_n_all","pt_n",100,0,30);
    TH1D* h_ee_eta_p_all = new TH1D("h_ee_eta_p_all","eta_p",100,-4,4);
    TH1D* h_ee_eta_n_all = new TH1D("h_ee_eta_n_all","eta_n",100,-4,4);
    TH1D* h_ee_phi_p_all = new TH1D("h_ee_phi_p_all","phi_p",100,-5,5);
    TH1D* h_ee_phi_n_all = new TH1D("h_ee_phi_n_all","phi_n",100,-5,5);
    TH1D* h_aco_all = new TH1D("h_aco_all","Acoplanarity",100,0,0.01);
    TH1D* h_pt_ee_all = new TH1D("h_pt_ee","Di-electron transverse momentum",100,0,5);
    TH1D* h_yy_all = new TH1D("h_yy_all","yy",100,-5,5);
    double lims[5] = {5,10,20,40,80};
    TH1D* h_xn0n08 = new TH1D("h_xn0n08","",4,lims);
    TH1D* h_0n0n08 = new TH1D("h_0n0n08","",4,lims);
    TH1D* h_xnxn08 = new TH1D("h_xnxn08","",4,lims);
    TH1D* h_xn0n16 = new TH1D("h_xn0n16","",4,lims);
    TH1D* h_0n0n16 = new TH1D("h_0n0n16","",4,lims);
    TH1D* h_xnxn16 = new TH1D("h_xnxn16","",4,lims);
    TH1D* h_xn0n24 = new TH1D("h_xn0n24","",4,lims);
    TH1D* h_xnxn24 = new TH1D("h_xnxn24","",4,lims);
    TH1D* h_0n0n24 = new TH1D("h_0n0n24","",4,lims);
    
    TH1D* h_xn0n08_c = new TH1D("h_xn0n08_c","",4,lims);
    TH1D* h_0n0n08_c = new TH1D("h_0n0n08_c","",4,lims);
    TH1D* h_xnxn08_c = new TH1D("h_xnxn08_c","",4,lims);
    TH1D* h_xn0n16_c = new TH1D("h_xn0n16_c","",4,lims);
    TH1D* h_0n0n16_c = new TH1D("h_0n0n16_c","",4,lims);
    TH1D* h_xnxn16_c = new TH1D("h_xnxn16_c","",4,lims);
    TH1D* h_xn0n24_c = new TH1D("h_xn0n24_c","",4,lims);
    TH1D* h_xnxn24_c = new TH1D("h_xnxn24_c","",4,lims);
    TH1D* h_0n0n24_c = new TH1D("h_0n0n24_c","",4,lims);
    
    TH1D* h_zdc_check_08 = new TH1D("h_zdc_check_08","",4,lims);
    TH1D* h_zdc_check_16 = new TH1D("h_zdc_check_16","",4,lims);
    TH1D* h_zdc_check_24 = new TH1D("h_zdc_check_24","",4,lims);
    //rapidity
    //rapidity_error
   	
   	
   	//rapidity: yll_ee
   	
   	
   	//open files
   	/* 
    TFile *mc_file_high = new TFile("output8M.root");
    TFile *mc_file_low = new TFile("outputP6M8.root");
    TFile *file_data = new TFile("outputdata.root");
    */
    /*
    TFile *mc_file_high = new TFile("output15M.root");
    TFile *mc_file_low = new TFile("outputP6M15.root");
    TFile *file_data = new TFile("outputdata.root");
    */
    TFile *mc_file_high = new TFile("output15Mv1.root");
    TFile *mc_file_low = new TFile("output7M15.root");
    TFile *mc_file_add = new TFile("output4p5M7.root");
    TFile *file_data = new TFile("outputdata.root");
    
    TFile *bkg = new TFile("fractionsStat.root");
    TFile *bkg_tight = new TFile("fractionsStatTight.root");
    
    //get histograms from files
    
    TH1D *h_ee_pt_p_high = (TH1D*)mc_file_high->Get("h_pt_p");
    TH1D *h_ee_pt_p_low = (TH1D*)mc_file_low->Get("h_pt_p");
    TH1D *h_ee_pt_p_dat = (TH1D*)file_data->Get("h_pt_p");
    
    TH1D *h_ee_pt_n_high = (TH1D*)mc_file_high->Get("h_pt_n");
    TH1D *h_ee_pt_n_low = (TH1D*)mc_file_low->Get("h_pt_n");
    TH1D *h_ee_pt_n_dat = (TH1D*)file_data->Get("h_pt_n");
    
    TH1D *h_ee_invm_high = (TH1D*)mc_file_high->Get("h_invm");
    TH1D *h_ee_invm_low = (TH1D*)mc_file_low->Get("h_invm");
    TH1D *h_ee_invm_dat = (TH1D*)file_data->Get("h_invm");
    
    TH1D *h_ee_eta_p_high = (TH1D*)mc_file_high->Get("h_eta_p");
    TH1D *h_ee_eta_p_low = (TH1D*)mc_file_low->Get("h_eta_p");
    TH1D *h_ee_eta_p_dat = (TH1D*)file_data->Get("h_eta_p");
    
    TH1D *h_ee_eta_n_high = (TH1D*)mc_file_high->Get("h_eta_n");
    TH1D *h_ee_eta_n_low = (TH1D*)mc_file_low->Get("h_eta_n");
    TH1D *h_ee_eta_n_dat = (TH1D*)file_data->Get("h_eta_n");
    
    TH1D *h_ee_phi_p_high = (TH1D*)mc_file_high->Get("h_phi_p");
    TH1D *h_ee_phi_p_low = (TH1D*)mc_file_low->Get("h_phi_p");
    TH1D *h_ee_phi_p_dat = (TH1D*)file_data->Get("h_phi_p");
    
    TH1D *h_ee_phi_n_high = (TH1D*)mc_file_high->Get("h_phi_n");
    TH1D *h_ee_phi_n_low = (TH1D*)mc_file_low->Get("h_phi_n");
    TH1D *h_ee_phi_n_dat = (TH1D*)file_data->Get("h_phi_n");
    
    TH1D *h_aco_high = (TH1D*)mc_file_high->Get("h_aco");
    TH1D *h_aco_low = (TH1D*)mc_file_low->Get("h_aco");
    TH1D *h_aco_dat = (TH1D*)file_data->Get("h_aco");
    
    TH1D *h_pt_ee_high = (TH1D*)mc_file_high->Get("Pt ee");
    TH1D *h_pt_ee_low = (TH1D*)mc_file_low->Get("Pt ee");
    TH1D *h_pt_ee_dat = (TH1D*)file_data->Get("Pt ee");
    
    TH1D *h_yy_high = (TH1D*)mc_file_high->Get("yy_ovf");
    TH1D *h_yy_low = (TH1D*)mc_file_low->Get("yy_ovf");
    TH1D *h_yy_dat = (TH1D*)file_data->Get("yy_ovf");
    
    TH2D* hh_zdc_ee = (TH2D*)file_data->Get("ZDC");
    TH2D* hh_pt_data = (TH2D*)file_data->Get("hh_pt1_pt2");
   	TH2D* hh_pt_low = (TH2D*)mc_file_low->Get("hh_pt1_pt2"); 
   	TH2D* hh_pt_high = (TH2D*)mc_file_high->Get("hh_pt1_pt2");
   	TH2D* hh_pt_mid = (TH2D*)mc_file_add->Get("hh_pt1_pt2");
    
    TH1D *h_m_ee_all08 = (TH1D*)file_data->Get("h_m_ee_all08");
    TH1D *h_m_ee_0n0n08 = (TH1D*)file_data->Get("h_m_ee_0n0n08");
    TH1D *h_m_ee_xn0n08 = (TH1D*)file_data->Get("h_m_ee_xn0n08");
    TH1D *h_m_ee_xnxn08 = (TH1D*)file_data->Get("h_m_ee_xnxn08");
    
    TH1D *h_m_ee_all16 = (TH1D*)file_data->Get("h_m_ee_all16");
    TH1D *h_m_ee_0n0n16 = (TH1D*)file_data->Get("h_m_ee_0n0n16");
    TH1D *h_m_ee_xn0n16= (TH1D*)file_data->Get("h_m_ee_xn0n16");
    TH1D *h_m_ee_xnxn16 = (TH1D*)file_data->Get("h_m_ee_xnxn16");
    
    TH1D *h_m_ee_all24 = (TH1D*)file_data->Get("h_m_ee_all24");
    TH1D *h_m_ee_0n0n24 = (TH1D*)file_data->Get("h_m_ee_0n0n24");
    TH1D *h_m_ee_xn0n24 = (TH1D*)file_data->Get("h_m_ee_xn0n24");
    TH1D *h_m_ee_xnxn24 = (TH1D*)file_data->Get("h_m_ee_xnxn24");
    
    TH1D *h_f0n0n = (TH1D*)file_data->Get("h_f0n0n");
    TH1D *h_fxn0n = (TH1D*)file_data->Get("h_fxn0n");
    TH1D *h_fxnxn = (TH1D*)file_data->Get("h_fxnxn");
    TH1D *h_f = (TH1D*)file_data->Get("h_f");
    
    TH1D *h_f0n0n_c = (TH1D*)h_f0n0n->Clone("h_f0n0n_c");
    TH1D *h_fxn0n_c = (TH1D*)h_fxn0n->Clone("h_fxn0n_c");
    TH1D *h_fxnxn_c = (TH1D*)h_fxnxn->Clone("h_fxnxn_c");
    
   	TH1D *h_f0n0n_o = (TH1D*)h_f0n0n->Clone("h_f0n0n_c");
    TH1D *h_fxn0n_o = (TH1D*)h_fxn0n->Clone("h_fxn0n_c");
    TH1D *h_fxnxn_o = (TH1D*)h_fxnxn->Clone("h_fxnxn_c");
    
    TH1D *h_ee_pt_p_add = (TH1D*)mc_file_add->Get("h_pt_p");  
    TH1D *h_ee_pt_n_add = (TH1D*)mc_file_add->Get("h_pt_n");
    TH1D *h_ee_invm_add = (TH1D*)mc_file_add->Get("h_invm");
    TH1D *h_ee_eta_p_add = (TH1D*)mc_file_add->Get("h_eta_p");
    TH1D *h_ee_eta_n_add = (TH1D*)mc_file_add->Get("h_eta_n");
    TH1D *h_ee_phi_p_add = (TH1D*)mc_file_add->Get("h_phi_p");
    TH1D *h_ee_phi_n_add = (TH1D*)mc_file_add->Get("h_phi_n");
    TH1D *h_aco_add = (TH1D*)mc_file_add->Get("h_aco"); 
    TH1D *h_pt_ee_add = (TH1D*)mc_file_add->Get("Pt ee");
    TH1D *h_yy_add = (TH1D*)mc_file_add->Get("yy_ovf");

    
    
    TH2D *fractions0n0n = (TH2D*)bkg->Get("allFractions0n0n_0");
    TH2D *fractionsxn0n = (TH2D*)bkg->Get("allFractionsXn0n_0");
    TH2D *fractionsxnxn = (TH2D*)bkg->Get("allFractionsXnXn_0");
    TH2D *fractionsinc = (TH2D*)bkg->Get("allFractionsInclusive_0");
    
    TH2D *Tfractions0n0n = (TH2D*)bkg_tight->Get("allFractions0n0n_0");
    TH2D *Tfractionsxn0n = (TH2D*)bkg_tight->Get("allFractionsXn0n_0");
    TH2D *Tfractionsxnxn = (TH2D*)bkg_tight->Get("allFractionsXnXn_0");
    TH2D *Tfractionsinc = (TH2D*)bkg_tight->Get("allFractionsInclusive_0");
    
    
    
    gStyle->SetErrorX(0);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickY(1);
    gStyle->SetPadTickX(1);
    
    //calculate weights
           
    double ZDC_factor = 1;   //0.7, 0.642
    /*
    double high_weight = 116*1440*ZDC_factor/500000.0; //N=100000 sig15=18.6 
    double low_weight = 419*1440*ZDC_factor/1500000.0; //sig3p6M15=518.3
    */
    /*
    double high_weight = 18.6*1440*ZDC_factor/100000.0; //N=100000 sig15=18.6 
    double low_weight = 518.3*1440*ZDC_factor/100000.0; //sig3p6M15=518.3
    */
    
    double high_weight = 17.83*1440*ZDC_factor/1000000.0; //N=1000000 sig15=17.83 
    double low_weight = 121.77*1440*ZDC_factor/4000000.0; //N=4000000 sig=121.77
    double add_weight = 156.55*1440*ZDC_factor/4000000.0; //N=4000000 sig=156.55
    
    cout << high_weight << " " << low_weight << " " << add_weight << endl;
    //scale histograms according to weights
    
    h_ee_pt_p_high->Scale(high_weight);
    h_ee_pt_p_low->Scale(low_weight);
    h_ee_invm_high->Scale(high_weight);
    h_ee_invm_low->Scale(low_weight);
    h_ee_pt_n_high->Scale(high_weight);
    h_ee_pt_n_low->Scale(low_weight);
    h_ee_eta_p_high->Scale(high_weight);
    h_ee_eta_p_low->Scale(low_weight);
    h_ee_eta_n_high->Scale(high_weight);
    h_ee_eta_n_low->Scale(low_weight);
    h_ee_phi_p_high->Scale(high_weight);
    h_ee_phi_p_low->Scale(low_weight);
    h_ee_phi_n_high->Scale(high_weight);
    h_ee_phi_n_low->Scale(low_weight);
    h_aco_high->Scale(high_weight);
    h_aco_low->Scale(low_weight);
    h_pt_ee_high->Scale(high_weight);
    h_pt_ee_low->Scale(low_weight);
    h_yy_high->Scale(high_weight);
    h_yy_low->Scale(low_weight);
    
    h_ee_pt_p_add->Scale(add_weight);
    h_ee_invm_add->Scale(add_weight);
    h_ee_pt_n_add->Scale(add_weight);
    h_ee_eta_p_add->Scale(add_weight);
    h_ee_eta_n_add->Scale(add_weight);
    h_ee_phi_p_add->Scale(add_weight);
    h_ee_phi_n_add->Scale(add_weight);
    h_aco_add->Scale(add_weight);
    h_pt_ee_add->Scale(add_weight);
    h_yy_add->Scale(add_weight);
    
    //divide rapidity histograms
    
    for (int i = 1 ;i<5;++i){
    	h_m_ee_0n0n08->SetBinContent(i,h_m_ee_0n0n08->GetBinContent(i)*(1-fractions0n0n->GetBinContent(i,1)));
    	h_m_ee_0n0n16->SetBinContent(i,h_m_ee_0n0n16->GetBinContent(i)*(1-fractions0n0n->GetBinContent(i,2)));
    	h_m_ee_0n0n24->SetBinContent(i,h_m_ee_0n0n24->GetBinContent(i)*(1-fractions0n0n->GetBinContent(i,3)));
    	
    	h_m_ee_xn0n08->SetBinContent(i,h_m_ee_xn0n08->GetBinContent(i)*(1-fractionsxn0n->GetBinContent(i,1)));
    	h_m_ee_xn0n16->SetBinContent(i,h_m_ee_xn0n16->GetBinContent(i)*(1-fractionsxn0n->GetBinContent(i,2)));
    	h_m_ee_xn0n24->SetBinContent(i,h_m_ee_xn0n24->GetBinContent(i)*(1-fractionsxn0n->GetBinContent(i,3)));
    	
    	h_m_ee_xnxn08->SetBinContent(i,h_m_ee_xnxn08->GetBinContent(i)*(1-fractionsxnxn->GetBinContent(i,1)));
    	h_m_ee_xnxn16->SetBinContent(i,h_m_ee_xnxn16->GetBinContent(i)*(1-fractionsxnxn->GetBinContent(i,2)));
    	h_m_ee_xnxn24->SetBinContent(i,h_m_ee_xnxn24->GetBinContent(i)*(1-fractionsxnxn->GetBinContent(i,3)));
    	
    	h_m_ee_all08->SetBinContent(i,h_m_ee_all08->GetBinContent(i)*(1-fractionsinc->GetBinContent(i,1)));
    	h_m_ee_all16->SetBinContent(i,h_m_ee_all16->GetBinContent(i)*(1-fractionsinc->GetBinContent(i,2)));
    	h_m_ee_all24->SetBinContent(i,h_m_ee_all24->GetBinContent(i)*(1-fractionsinc->GetBinContent(i,3)));
    
    }
    
    
    
    h_xn0n08->Divide(h_m_ee_xn0n08,h_m_ee_all08,1,1,"b");
    h_xnxn08->Divide(h_m_ee_xnxn08,h_m_ee_all08,1,1,"b");
    h_0n0n08->Divide(h_m_ee_0n0n08,h_m_ee_all08,1,1,"b");
    h_xn0n16->Divide(h_m_ee_xn0n16,h_m_ee_all16,1,1,"b");
    h_xnxn16->Divide(h_m_ee_xnxn16,h_m_ee_all16,1,1,"b");
    h_0n0n16->Divide(h_m_ee_0n0n16,h_m_ee_all16,1,1,"b");
    h_xn0n24->Divide(h_m_ee_xn0n24,h_m_ee_all24,1,1,"b");
    h_xnxn24->Divide(h_m_ee_xnxn24,h_m_ee_all24,1,1,"b");
    h_0n0n24->Divide(h_m_ee_0n0n24,h_m_ee_all24,1,1,"b");
    
    

    
    
    //odjąć tło <<<<<<<<<-------------------------------------------------
    //dysocjacja 
    


    h_f0n0n->SetBinContent(1,h_f0n0n->GetBinContent(1)*(1-fractions0n0n->GetBinContent(5,4)));
    h_fxn0n->SetBinContent(1,h_fxn0n->GetBinContent(1)*(1-fractionsxn0n->GetBinContent(5,4)));
    h_fxnxn->SetBinContent(1,h_fxnxn->GetBinContent(1)*(1-fractionsxnxn->GetBinContent(5,4)));
    h_f->SetBinContent(1,h_f->GetBinContent(1)*(1-fractionsinc->GetBinContent(5,4)));
    
    h_f0n0n->SetBinContent(2,h_f0n0n->GetBinContent(2)*(1-Tfractions0n0n->GetBinContent(4,4)));
    h_fxn0n->SetBinContent(2,h_fxn0n->GetBinContent(2)*(1-Tfractionsxn0n->GetBinContent(4,4)));
    h_fxnxn->SetBinContent(2,h_fxnxn->GetBinContent(2)*(1-Tfractionsxnxn->GetBinContent(4,4)));
    h_f->SetBinContent(2,h_f->GetBinContent(2)*(1-Tfractionsinc->GetBinContent(4,4)));
    

    //h_f0n0n *(1-liczba)
    //h_f *(1-inna_liczba)
    
    
    h_f0n0n->Divide(h_f0n0n,h_f,1,1,"b");
    h_fxn0n->Divide(h_fxn0n,h_f,1,1,"b");
    h_fxnxn->Divide(h_fxnxn,h_f,1,1,"b");

    

    
    //h_m_ee_xnxn->Divide(h_m_ee_all);
    //h_m_ee_xn0n->Divide(h_m_ee_all);
    
    
    //TAsymmError
    
    //pileup corrections
    double f0n0n08, f0n0n16, f0n0n24, fxn0n08, fxn0n16, fxn0n24, fxnxn08, fxnxn16, fxnxn24;
    double a00, a10, a11, a20, a21, a22;
    
    double e00, e01, e02, e10, e11, e12, e13, e20, e21, e22, e23, e24, error0, error1, error2;
    
    double ps= 6.362e-2, pm=1.978e-3; // err ps +0.38 -0.33, err pm +-0.011 2015: ps= 5.59e-2, pm=0.174e-2, ps = 6.64e-2, pm = 1.98e-3       ps =6.362e-2
    double ps_en = 0.367e-2, ps_ep = 0.429e-2, pm_e = 0.136e-3; //ps_en = 0.38e-2, ps_ep = 0.33e-2, pm_e = 0.011e-2;
    
    double f00, f0n, fnn;
    	
    a00 = -1/((ps-1)*(ps-1)*(pm-1));
    a10 = (ps*(ps*(pm-2)-2*pm+2))/((ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1));
    a11 = 1/((ps-1)*(pm-1));
    a20 = (ps*ps*ps*(pm-1)*(pm-1)+ps*ps*(-3*pm*pm+5*pm-1)+3*ps*(pm-1)*pm-(pm-1)*pm)/((ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1));
    a21 = (ps*(pm-1)-pm)/((pm-1)*(ps-1));
    a22 = 1;
    
   
    TGraphAsymmErrors* err_0n0n08 = new TGraphAsymmErrors(f0n0n08);
    TGraphAsymmErrors* err_xn0n08 = new TGraphAsymmErrors(fxn0n08);
    TGraphAsymmErrors* err_xnxn08 = new TGraphAsymmErrors(fxnxn08);
    TGraphAsymmErrors* err_0n0n16 = new TGraphAsymmErrors(f0n0n16);
    TGraphAsymmErrors* err_xn0n16 = new TGraphAsymmErrors(fxn0n16);
    TGraphAsymmErrors* err_xnxn16 = new TGraphAsymmErrors(fxnxn16);
    TGraphAsymmErrors* err_0n0n24 = new TGraphAsymmErrors(f0n0n24);
    TGraphAsymmErrors* err_xn0n24 = new TGraphAsymmErrors(fxn0n24);
    TGraphAsymmErrors* err_xnxn24 = new TGraphAsymmErrors(fxnxn24);
    
    TGraphAsymmErrors* o_0n0n08 = new TGraphAsymmErrors(f0n0n08);
    TGraphAsymmErrors* o_xn0n08 = new TGraphAsymmErrors(fxn0n08);
    TGraphAsymmErrors* o_xnxn08 = new TGraphAsymmErrors(fxnxn08);
    TGraphAsymmErrors* o_0n0n16 = new TGraphAsymmErrors(f0n0n16);
    TGraphAsymmErrors* o_xn0n16 = new TGraphAsymmErrors(fxn0n16);
    TGraphAsymmErrors* o_xnxn16 = new TGraphAsymmErrors(fxnxn16);
    TGraphAsymmErrors* o_0n0n24 = new TGraphAsymmErrors(f0n0n24);
    TGraphAsymmErrors* o_xn0n24 = new TGraphAsymmErrors(fxn0n24);
    TGraphAsymmErrors* o_xnxn24 = new TGraphAsymmErrors(fxnxn24);
    
    
    double err_0n0n08_x[4], err_0n0n08_y[4],err_0n0n08_eyl[4],err_0n0n08_eyh[4],err_x0[4]={0};
    double err_xn0n08_x[4], err_xn0n08_y[4],err_xn0n08_eyl[4],err_xn0n08_eyh[4];
	double err_xnxn08_x[4], err_xnxn08_y[4],err_xnxn08_eyl[4],err_xnxn08_eyh[4]; 
	double err_0n0n16_x[4], err_0n0n16_y[4],err_0n0n16_eyl[4],err_0n0n16_eyh[4];
    double err_xn0n16_x[4], err_xn0n16_y[4],err_xn0n16_eyl[4],err_xn0n16_eyh[4];
	double err_xnxn16_x[4], err_xnxn16_y[4],err_xnxn16_eyl[4],err_xnxn16_eyh[4];
	double err_0n0n24_x[4], err_0n0n24_y[4],err_0n0n24_eyl[4],err_0n0n24_eyh[4];
    double err_xn0n24_x[4], err_xn0n24_y[4],err_xn0n24_eyl[4],err_xn0n24_eyh[4];
	double err_xnxn24_x[4], err_xnxn24_y[4],err_xnxn24_eyl[4],err_xnxn24_eyh[4];     
    
    double o_0n0n08_x[4], o_0n0n08_y[4],o_0n0n08_eyl[4],o_0n0n08_eyh[4],o_x0[4]={0};
    double o_xn0n08_x[4], o_xn0n08_y[4],o_xn0n08_eyl[4],o_xn0n08_eyh[4];
	double o_xnxn08_x[4], o_xnxn08_y[4],o_xnxn08_eyl[4],o_xnxn08_eyh[4]; 
	double o_0n0n16_x[4], o_0n0n16_y[4],o_0n0n16_eyl[4],o_0n0n16_eyh[4];
    double o_xn0n16_x[4], o_xn0n16_y[4],o_xn0n16_eyl[4],o_xn0n16_eyh[4];
	double o_xnxn16_x[4], o_xnxn16_y[4],o_xnxn16_eyl[4],o_xnxn16_eyh[4];
	double o_0n0n24_x[4], o_0n0n24_y[4],o_0n0n24_eyl[4],o_0n0n24_eyh[4];
    double o_xn0n24_x[4], o_xn0n24_y[4],o_xn0n24_eyl[4],o_xn0n24_eyh[4];
	double o_xnxn24_x[4], o_xnxn24_y[4],o_xnxn24_eyl[4],o_xnxn24_eyh[4]; 
    
    
    
    
   	double error_x;
   	err_0n0n08_x[0]=7.5;
   	err_0n0n08_x[1]=15;
   	err_0n0n08_x[2]=30;
   	err_0n0n08_x[3]=60;
   	
   	o_0n0n08_x[0]=7.5;
   	o_0n0n08_x[1]=15;
   	o_0n0n08_x[2]=30;
   	o_0n0n08_x[3]=60;
   	
    for(int i=1;i<=4;++i){
    
    	
    	
    
    	f0n0n08 = h_0n0n08->GetBinContent(i);
    	f0n0n16 = h_0n0n16->GetBinContent(i);
    	f0n0n24 = h_0n0n24->GetBinContent(i);
    	fxn0n08 = h_xn0n08->GetBinContent(i);
    	fxn0n16 = h_xn0n16->GetBinContent(i);
		fxn0n24 = h_xn0n24->GetBinContent(i);
		fxnxn08 = h_xnxn08->GetBinContent(i);
    	fxnxn16 = h_xnxn16->GetBinContent(i);
    	fxnxn24 = h_xnxn24->GetBinContent(i);
    	
    	o_0n0n08_y[i-1]=f0n0n08;
    	o_xn0n08_y[i-1]=fxn0n08;
    	o_xnxn08_y[i-1]=fxnxn08;
    	o_0n0n16_y[i-1]=f0n0n16;
    	o_xn0n16_y[i-1]=fxn0n16;
    	o_xnxn16_y[i-1]=fxnxn16;
    	o_0n0n24_y[i-1]=f0n0n24;
    	o_xn0n24_y[i-1]=fxn0n24;
    	o_xnxn24_y[i-1]=fxnxn24;

    	h_0n0n08_c->SetBinContent(i,f0n0n08*a00);
    	h_0n0n16_c->SetBinContent(i,f0n0n16*a00);   	
    	h_0n0n24_c->SetBinContent(i,f0n0n24*a00);
    	
    	h_xn0n08_c->SetBinContent(i,a10*f0n0n08+a11*fxn0n08);
    	h_xn0n16_c->SetBinContent(i,a10*f0n0n16+a11*fxn0n16);
    	h_xn0n24_c->SetBinContent(i,a10*f0n0n24+a11*fxn0n24);
    	
    	h_xnxn08_c->SetBinContent(i,a20*f0n0n08+a21*fxn0n08+a22*fxnxn08);
    	h_xnxn16_c->SetBinContent(i,a20*f0n0n16+a21*fxn0n16+a22*fxnxn16);  	
    	h_xnxn24_c->SetBinContent(i,a20*f0n0n24+a21*fxn0n24+a22*fxnxn24);

    	/*
    	h_0n0n16_c->GetBinContent(i);
    	h_0n0n24_c->GetBinContent(i);
    	h_xn0n08_c->GetBinContent(i);
    	h_xn0n16_c->GetBinContent(i);
    	h_xn0n24_c->GetBinContent(i);
    	h_xnxn08_c->GetBinContent(i);
    	h_xnxn16_c->GetBinContent(i);  	
    	h_xnxn24_c->GetBinContent(i);
    	*/
    	
    	//calculate errors 08
    	f00 = h_0n0n08->GetBinError(i);
		f0n = h_xn0n08->GetBinError(i);
		fnn = h_xnxn08->GetBinError(i);
		
		o_0n0n08_eyl[i-1]=f00;
		o_xn0n08_eyl[i-1]=f0n;
		o_xnxn08_eyl[i-1]=fnn;
		o_0n0n08_eyh[i-1]=f00;
		o_xn0n08_eyh[i-1]=f0n;
		o_xnxn08_eyh[i-1]=fnn;

    	e00 = 2*f00/((ps-1)*(ps-1)*(ps-1)*(pm-1));
    	e01 = 1*f00/((ps-1)*(ps-1)*(pm-1)*(pm-1));
    	e02 = a00;
    
    	e10 = (f00*(ps*ps*(-(pm-2))+2*ps*pm+2*(pm-1))-f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1));
    	e11 = (f00*ps*(ps*(-pm)+3*ps+2*pm-2)-f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1)*(pm-1));
    	e12 = a10;
    	e13 = a11;
    
    	e20 = (f00*ps*(ps*(pm-2)-4*pm+2)+f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1));
    	e21 = (f00*(ps*ps*(pm-3)-3*ps*(pm-1)+pm+1)+f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1)*(pm-1));
    	e22 = a20;
    	e23 = a21;
    	e24 = a22;
    	
    	err_0n0n08_y[i-1]=h_0n0n08_c->GetBinContent(i);
    	err_0n0n16_y[i-1]=h_0n0n16_c->GetBinContent(i);
		err_0n0n24_y[i-1]=h_0n0n24_c->GetBinContent(i);
    	err_xn0n08_y[i-1]=h_xn0n08_c->GetBinContent(i);
    	err_xn0n16_y[i-1]=h_xn0n16_c->GetBinContent(i);
		err_xn0n24_y[i-1]=h_xn0n24_c->GetBinContent(i);
		err_xnxn08_y[i-1]=h_xnxn08_c->GetBinContent(i); 
    	err_xnxn16_y[i-1]=h_xnxn16_c->GetBinContent(i);
    	err_xnxn24_y[i-1]=h_xnxn24_c->GetBinContent(i); 
    	
   		
   		error_x=sqrt((e00*ps_ep)*(e00*ps_ep)+(e01*pm_e)*(e01*pm_e)+(e02*f00)*(e02*f00));
   		err_0n0n08->SetPointEYhigh(i,error_x);
   		err_0n0n08_eyh[i-1]=error_x;
   		
    	error_x=sqrt((e00*ps_en)*(e00*ps_en)+(e01*pm_e)*(e01*pm_e)+(e02*f00)*(e02*f00));
    	err_0n0n08->SetPointEYlow(i,error_x);
    	err_0n0n08_eyl[i-1]=error_x;
    	
    	
    	
    	error_x=sqrt((e10*ps_ep)*(e10*ps_ep)+(e11*pm_e)*(e11*pm_e)+(e12*f00)*(e12*f00)+(e13*f0n)*(e13*f0n));
    	err_xn0n08->SetPointEYhigh(i,error_x);
    	err_xn0n08_eyh[i-1]=error_x;
    	error_x=sqrt((e10*ps_en)*(e10*ps_en)+(e11*pm_e)*(e11*pm_e)+(e12*f00)*(e12*f00)+(e13*f0n)*(e13*f0n));
    	err_xn0n08->SetPointEYlow(i,error_x);
    	err_xn0n08_eyl[i-1]=error_x;
    	
    	error_x=sqrt((e20*ps_ep)*(e20*ps_ep)+(e21*pm_e)*(e21*pm_e)+(e22*f00)*(e22*f00)+(e23*f0n)*(e23*f0n)+(e24*fnn)*(e24*fnn));
    	err_xnxn08->SetPointEYhigh(i,error_x);
    	err_xnxn08_eyh[i-1]=error_x;
    	error_x=sqrt((e20*ps_en)*(e20*ps_en)+(e21*pm_e)*(e21*pm_e)+(e22*f00)*(e22*f00)+(e23*f0n)*(e23*f0n)+(e24*fnn)*(e24*fnn));
    	err_xnxn08->SetPointEYlow(i,error_x);
    	err_xnxn08_eyl[i-1]=error_x;
    	
    	//errors 16
    	
    	f00 = h_0n0n16->GetBinError(i);
		f0n = h_xn0n16->GetBinError(i);
		fnn = h_xnxn16->GetBinError(i);
		
		o_0n0n16_eyl[i-1]=f00;
		o_xn0n16_eyl[i-1]=f0n;
		o_xnxn16_eyl[i-1]=fnn;
		o_0n0n16_eyh[i-1]=f00;
		o_xn0n16_eyh[i-1]=f0n;
		o_xnxn16_eyh[i-1]=fnn;
    	
    	e00 = 2*f00/((ps-1)*(ps-1)*(ps-1)*(pm-1));
    	e01 = 1*f00/((ps-1)*(ps-1)*(pm-1)*(pm-1));
    	e02 = a00;
    
    	e10 = (f00*(ps*ps*(-(pm-2))+2*ps*pm+2*(pm-1))-f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1));
    	e11 = (f00*ps*(ps*(-pm)+3*ps+2*pm-2)-f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1)*(pm-1));
    	e12 = a10;
    	e13 = a11;
    
    	e20 = (f00*ps*(ps*(pm-2)-4*pm+2)+f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1));
    	e21 = (f00*(ps*ps*(pm-3)-3*ps*(pm-1)+pm+1)+f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1)*(pm-1));
    	e22 = a20;
    	e23 = a21;
    	e24 = a22;
   		
   		error_x=sqrt((e00*ps_ep)*(e00*ps_ep)+(e01*pm_e)*(e01*pm_e)+(e02*f00)*(e02*f00));
   		err_0n0n16->SetPointEYhigh(i,error_x);
   		err_0n0n16_eyh[i-1]=error_x;
    	error_x=sqrt((e00*ps_en)*(e00*ps_en)+(e01*pm_e)*(e01*pm_e)+(e02*f00)*(e02*f00));
    	err_0n0n16->SetPointEYlow(i,error_x);
    	err_0n0n16_eyl[i-1]=error_x;
    	
    	error_x=sqrt((e10*ps_ep)*(e10*ps_ep)+(e11*pm_e)*(e11*pm_e)+(e12*f00)*(e12*f00)+(e13*f0n)*(e13*f0n));
    	err_xn0n16->SetPointEYhigh(i,error_x);
    	err_xn0n16_eyh[i-1]=error_x;
    	error_x=sqrt((e10*ps_en)*(e10*ps_en)+(e11*pm_e)*(e11*pm_e)+(e12*f00)*(e12*f00)+(e13*f0n)*(e13*f0n));
    	err_xn0n16->SetPointEYlow(i,error_x);
    	err_xn0n16_eyl[i-1]=error_x;
    	
    	error_x=sqrt((e20*ps_ep)*(e20*ps_ep)+(e21*pm_e)*(e21*pm_e)+(e22*f00)*(e22*f00)+(e23*f0n)*(e23*f0n)+(e24*fnn)*(e24*fnn));
    	err_xnxn16->SetPointEYhigh(i,error_x);
    	err_xnxn16_eyh[i-1]=error_x;
    	error_x=sqrt((e20*ps_en)*(e20*ps_en)+(e21*pm_e)*(e21*pm_e)+(e22*f00)*(e22*f00)+(e23*f0n)*(e23*f0n)+(e24*fnn)*(e24*fnn));
    	err_xnxn16->SetPointEYlow(i,error_x);
    	err_xnxn16_eyl[i-1]=error_x;
    	
    	//errors 24
    	
    	f00 = h_0n0n24->GetBinError(i);
		f0n = h_xn0n24->GetBinError(i);
		fnn = h_xnxn24->GetBinError(i);
		
		o_0n0n24_eyl[i-1]=f00;
		o_xn0n24_eyl[i-1]=f0n;
		o_xnxn24_eyl[i-1]=fnn;
		o_0n0n24_eyh[i-1]=f00;
		o_xn0n24_eyh[i-1]=f0n;
		o_xnxn24_eyh[i-1]=fnn;

    	e00 = 2*f00/((ps-1)*(ps-1)*(ps-1)*(pm-1));
    	e01 = 1*f00/((ps-1)*(ps-1)*(pm-1)*(pm-1));
    	e02 = a00;
    
    	e10 = (f00*(ps*ps*(-(pm-2))+2*ps*pm+2*(pm-1))-f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1));
    	e11 = (f00*ps*(ps*(-pm)+3*ps+2*pm-2)-f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1)*(pm-1));
    	e12 = a10;
    	e13 = a11;
    
    	e20 = (f00*ps*(ps*(pm-2)-4*pm+2)+f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1));
    	e21 = (f00*(ps*ps*(pm-3)-3*ps*(pm-1)+pm+1)+f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1)*(pm-1));
    	e22 = a20;
    	e23 = a21;
    	e24 = a22;
   		
   		error_x=sqrt((e00*ps_ep)*(e00*ps_ep)+(e01*pm_e)*(e01*pm_e)+(e02*f00)*(e02*f00));
   		err_0n0n24->SetPointEYhigh(i,error_x);
   		err_0n0n24_eyh[i-1]=error_x;
    	error_x=sqrt((e00*ps_en)*(e00*ps_en)+(e01*pm_e)*(e01*pm_e)+(e02*f00)*(e02*f00));
    	err_0n0n24->SetPointEYlow(i,error_x);
    	err_0n0n24_eyl[i-1]=error_x;
    	
    	error_x=sqrt((e10*ps_ep)*(e10*ps_ep)+(e11*pm_e)*(e11*pm_e)+(e12*f00)*(e12*f00)+(e13*f0n)*(e13*f0n));
    	err_xn0n24->SetPointEYhigh(i,error_x);
    	err_xn0n24_eyh[i-1]=error_x;
    	error_x=sqrt((e10*ps_en)*(e10*ps_en)+(e11*pm_e)*(e11*pm_e)+(e12*f00)*(e12*f00)+(e13*f0n)*(e13*f0n));
    	err_xn0n24->SetPointEYlow(i,error_x);
    	err_xn0n24_eyl[i-1]=error_x;
    	
    	error_x=sqrt((e20*ps_ep)*(e20*ps_ep)+(e21*pm_e)*(e21*pm_e)+(e22*f00)*(e22*f00)+(e23*f0n)*(e23*f0n)+(e24*fnn)*(e24*fnn));
    	err_xnxn24->SetPointEYhigh(i,error_x);
    	err_xnxn24_eyh[i-1]=error_x;
    	error_x=sqrt((e20*ps_en)*(e20*ps_en)+(e21*pm_e)*(e21*pm_e)+(e22*f00)*(e22*f00)+(e23*f0n)*(e23*f0n)+(e24*fnn)*(e24*fnn));
    	err_xnxn24->SetPointEYlow(i,error_x);
    	err_xnxn24_eyl[i-1]=error_x;
    	
    	
    }


	
	
	TGraphAsymmErrors* n_0n0n08 = new TGraphAsymmErrors(4,err_0n0n08_x,err_0n0n08_y,err_x0,err_x0,err_0n0n08_eyl,err_0n0n08_eyh);
	TGraphAsymmErrors* n_xn0n08 = new TGraphAsymmErrors(4,err_0n0n08_x,err_xn0n08_y,err_x0,err_x0,err_xn0n08_eyl,err_xn0n08_eyh);
	TGraphAsymmErrors* n_xnxn08 = new TGraphAsymmErrors(4,err_0n0n08_x,err_xnxn08_y,err_x0,err_x0,err_xnxn08_eyl,err_xnxn08_eyh);
	TGraphAsymmErrors* n_0n0n16 = new TGraphAsymmErrors(4,err_0n0n08_x,err_0n0n16_y,err_x0,err_x0,err_0n0n16_eyl,err_0n0n16_eyh);
	TGraphAsymmErrors* n_xn0n16 = new TGraphAsymmErrors(4,err_0n0n08_x,err_xn0n16_y,err_x0,err_x0,err_xn0n16_eyl,err_xn0n16_eyh);
	TGraphAsymmErrors* n_xnxn16 = new TGraphAsymmErrors(4,err_0n0n08_x,err_xnxn16_y,err_x0,err_x0,err_xnxn16_eyl,err_xnxn16_eyh);
	TGraphAsymmErrors* n_0n0n24 = new TGraphAsymmErrors(4,err_0n0n08_x,err_0n0n24_y,err_x0,err_x0,err_0n0n24_eyl,err_0n0n24_eyh);
	TGraphAsymmErrors* n_xn0n24 = new TGraphAsymmErrors(4,err_0n0n08_x,err_xn0n24_y,err_x0,err_x0,err_xn0n24_eyl,err_xn0n24_eyh);
	TGraphAsymmErrors* n_xnxn24 = new TGraphAsymmErrors(4,err_0n0n08_x,err_xnxn24_y,err_x0,err_x0,err_xnxn24_eyl,err_xnxn24_eyh);
	

	TGraphAsymmErrors* o_0n0n08_h = new TGraphAsymmErrors(4,o_0n0n08_x,o_0n0n08_y,err_x0,err_x0,o_0n0n08_eyl,o_0n0n08_eyh);
	TGraphAsymmErrors* o_0n0n16_h = new TGraphAsymmErrors(4,o_0n0n08_x,o_0n0n16_y,err_x0,err_x0,o_0n0n16_eyl,o_0n0n16_eyh);
	TGraphAsymmErrors* o_0n0n24_h = new TGraphAsymmErrors(4,o_0n0n08_x,o_0n0n24_y,err_x0,err_x0,o_0n0n24_eyl,o_0n0n24_eyh);
	TGraphAsymmErrors* o_xn0n08_h = new TGraphAsymmErrors(4,o_0n0n08_x,o_xn0n08_y,err_x0,err_x0,o_xn0n08_eyl,o_xn0n08_eyh);
	TGraphAsymmErrors* o_xn0n16_h = new TGraphAsymmErrors(4,o_0n0n08_x,o_xn0n16_y,err_x0,err_x0,o_xn0n16_eyl,o_xn0n16_eyh);
	TGraphAsymmErrors* o_xn0n24_h = new TGraphAsymmErrors(4,o_0n0n08_x,o_xn0n24_y,err_x0,err_x0,o_xn0n24_eyl,o_xn0n24_eyh);
	TGraphAsymmErrors* o_xnxn08_h = new TGraphAsymmErrors(4,o_0n0n08_x,o_xnxn08_y,err_x0,err_x0,o_xnxn08_eyl,o_xnxn08_eyh);
	TGraphAsymmErrors* o_xnxn16_h = new TGraphAsymmErrors(4,o_0n0n08_x,o_xnxn16_y,err_x0,err_x0,o_xnxn16_eyl,o_xnxn16_eyh);
	TGraphAsymmErrors* o_xnxn24_h = new TGraphAsymmErrors(4,o_0n0n08_x,o_xnxn24_y,err_x0,err_x0,o_xnxn24_eyl,o_xnxn24_eyh);

 


	double of_0n0n_x[2], of_0n0n_y[2],of_0n0n_eyl[2],of_0n0n_eyh[2];
    double of_xn0n_x[2], of_xn0n_y[2],of_xn0n_eyl[2],of_xn0n_eyh[2];
    double of_xnxn_x[2], of_xnxn_y[2],of_xnxn_eyl[2],of_xnxn_eyh[2];
   
   	of_0n0n_x[0]=0.5;
    of_0n0n_x[1]=1.5;
	of_xn0n_x[0]=0.5;
    of_xn0n_x[1]=1.5;
	of_xnxn_x[0]=0.5;
    of_xnxn_x[1]=1.5;
   
    double c0n0n, cxn0n, cxnxn;
    
    h_f0n0n_o->Divide(h_f0n0n_c,h_f,1,1,"b");
    h_fxn0n_o->Divide(h_fxn0n_c,h_f,1,1,"b");
    h_fxnxn_o->Divide(h_fxnxn_c,h_f,1,1,"b");
    
    for(int i=1; i<=2;++i){
    	c0n0n = h_f0n0n_c->GetBinContent(i);
    	cxn0n = h_fxn0n_c->GetBinContent(i);
    	cxnxn = h_fxnxn_c->GetBinContent(i);
    	
    	of_0n0n_y[i-1]=h_f0n0n_o->GetBinContent(i);
    	of_xn0n_y[i-1]=h_fxn0n_o->GetBinContent(i);	
    	of_xnxn_y[i-1]=h_fxnxn_o->GetBinContent(i);	
    	
    	of_0n0n_eyl[i-1]=h_f0n0n_o->GetBinError(i);
    	of_xn0n_eyl[i-1]=h_fxn0n_o->GetBinError(i);
    	of_xnxn_eyl[i-1]=h_fxnxn_o->GetBinError(i);
    	
    	of_0n0n_eyh[i-1]=h_f0n0n_o->GetBinError(i);
    	of_xn0n_eyh[i-1]=h_fxn0n_o->GetBinError(i);
    	of_xnxn_eyh[i-1]=h_fxnxn_o->GetBinError(i);
    	
    	h_f0n0n_c->SetBinContent(i,c0n0n*a00);
    	h_fxn0n_c->SetBinContent(i,a10*c0n0n+a11*cxn0n);
    	h_fxnxn_c->SetBinContent(i,a20*c0n0n+a21*cxn0n+a22*cxnxn);
    	
    	
    	

    }
    cout << of_0n0n_y[0] << " " << of_0n0n_y[1] << endl;
    h_f0n0n_c->Divide(h_f0n0n_c,h_f,1,1,"b");
    h_fxn0n_c->Divide(h_fxn0n_c,h_f,1,1,"b");
    h_fxnxn_c->Divide(h_fxnxn_c,h_f,1,1,"b");
    cout << "pre" << endl;
    cout << "0n0n08" <<endl;
    for (int i =0;i<4;++i){
    	cout << o_0n0n08_y[i] << endl;//o_0n0n08_eyl[i] << " " << o_0n0n08_eyh[i] << endl;
    }
    cout << "xn0n08" << endl;
    for (int i =0;i<4;++i){
    	cout <<o_xn0n08_y[i] << endl;
    }
    
    cout << "xnxn08" << endl;
    for (int i =0;i<4;++i){
    	cout << o_xnxn08_y[i] << endl;
    }
    cout<< "0n0n16" << endl;
    for (int i =0;i<4;++i){
    	cout << o_0n0n16_y[i] << endl;
    }
    cout << "xn0n16" << endl;
    for (int i =0;i<4;++i){
    	cout << o_xn0n16_y[i]<< endl;
    }
    
    cout << "xnxn16" << endl;
    for (int i =0;i<4;++i){
    	cout << o_xnxn16_y[i] << endl;
    }
    cout << "0n0n24" << endl;
    for (int i =0;i<4;++i){
    	cout << o_0n0n24_y[i]<< endl;
    }
    cout << "xn0n24" << endl;
    for (int i =0;i<4;++i){
    	cout << o_xn0n24_y[i] << endl;
    }
    
    cout << "xnxn24" << endl;
    for (int i =0;i<4;++i){
    	cout << o_xnxn24_y[i]<< endl;
    }
        cout << "post" << endl;
    
    cout << "0n0n08" <<endl;
    for (int i =0;i<4;++i){
    	cout << err_0n0n08_y[i] << endl;// err_0n0n08_eyl[i] << " " << err_0n0n08_eyh[i] << endl;
    }
    cout << "xn0n08" << endl;
    for (int i =0;i<4;++i){
    	cout <<err_xn0n08_y[i] << endl;
    }
    

    
    cout << "xnxn08" << endl;
    for (int i =0;i<4;++i){
    	cout << err_xnxn08_y[i] << endl;
    }
    cout<< "0n0n16" << endl;
    for (int i =0;i<4;++i){
    	cout << err_0n0n16_y[i] << endl;
    }
    cout << "xn0n16" << endl;
    for (int i =0;i<4;++i){
    	cout << err_xn0n16_y[i]<< endl;
    }
    
    cout << "xnxn16" << endl;
    for (int i =0;i<4;++i){
    	cout << err_xnxn16_y[i] << endl;
    }
    cout << "0n0n24" << endl;
    for (int i =0;i<4;++i){
    	cout << err_0n0n24_y[i]<< endl;
    }
    cout << "xn0n24" << endl;
    for (int i =0;i<4;++i){
    	cout << err_xn0n24_y[i] << endl;
    }
    
    cout << "xnxn24" << endl;
    for (int i =0;i<4;++i){
    	cout << err_xnxn24_y[i]<< endl;
    }
    
    double f_0n0n_x[2], f_0n0n_y[2],f_0n0n_eyl[2],f_0n0n_eyh[2];
    double f_xn0n_x[2], f_xn0n_y[2],f_xn0n_eyl[2],f_xn0n_eyh[2];
    double f_xnxn_x[2], f_xnxn_y[2],f_xnxn_eyl[2],f_xnxn_eyh[2];
    
    f_0n0n_x[0]=0.5;
    f_0n0n_x[1]=1.5;
	f_xn0n_x[0]=0.5;
    f_xn0n_x[1]=1.5;
	f_xnxn_x[0]=0.5;
    f_xnxn_x[1]=1.5;

    TGraphAsymmErrors* err_0n0n_c = new TGraphAsymmErrors(h_f0n0n_c);
    TGraphAsymmErrors* err_xn0n_c = new TGraphAsymmErrors(h_fxn0n_c);
    TGraphAsymmErrors* err_xnxn_c = new TGraphAsymmErrors(h_fxnxn_c);
    
    cout << "Asymerror: " << h_f0n0n_c->GetBinError(2) << endl;
    
    for(int i=1; i<=2;++i){
    	f00 = h_f0n0n_c->GetBinError(i);
		f0n = h_fxn0n_c->GetBinError(i);
		fnn = h_fxnxn_c->GetBinError(i);
		
    	f_0n0n_y[i-1]=h_f0n0n_c->GetBinContent(i);
    	f_xn0n_y[i-1]=h_fxn0n_c->GetBinContent(i);
    	f_xnxn_y[i-1]=h_fxnxn_c->GetBinContent(i);
    	
    	e00 = 2*f00/((ps-1)*(ps-1)*(ps-1)*(pm-1));
    	e01 = 1*f00/((ps-1)*(ps-1)*(pm-1)*(pm-1));
    	e02 = a00;
    
    	e10 = (f00*(ps*ps*(-(pm-2))+2*ps*pm+2*(pm-1))-f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1));
    	e11 = (f00*ps*(ps*(-pm)+3*ps+2*pm-2)-f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1)*(pm-1));
    	e12 = a10;
    	e13 = a11;
    
    	e20 = (f00*ps*(ps*(pm-2)-4*pm+2)+f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1));
    	e21 = (f00*(ps*ps*(pm-3)-3*ps*(pm-1)+pm+1)+f0n*(ps-1)*(ps-1)*(pm-1))/((ps-1)*(ps-1)*(ps-1)*(pm-1)*(pm-1)*(pm-1));
    	e22 = a20;
    	e23 = a21;
    	e24 = a22;
   		
   		error_x=sqrt((e00*ps_ep)*(e00*ps_ep)+(e01*pm_e)*(e01*pm_e)+(e02*f00)*(e02*f00));
   		//cout << "error_x_high: " << error_x << endl;
   		err_0n0n_c->SetPointEYhigh(i,error_x);
   		f_0n0n_eyh[i-1]=error_x;
   		
    	error_x=sqrt((e00*ps_en)*(e00*ps_en)+(e01*pm_e)*(e01*pm_e)+(e02*f00)*(e02*f00));
    	err_0n0n_c->SetPointEYlow(i,error_x);
    	f_0n0n_eyl[i-1]=error_x;
    	//cout << "error_x_low: " << error_x << endl;
    	
    	error_x=sqrt((e10*ps_ep)*(e10*ps_ep)+(e11*pm_e)*(e11*pm_e)+(e12*f00)*(e12*f00)+(e13*f0n)*(e13*f0n));
    	err_xn0n_c->SetPointEYhigh(i,error_x);
    	f_xn0n_eyh[i-1]=error_x;
    	error_x=sqrt((e10*ps_en)*(e10*ps_en)+(e11*pm_e)*(e11*pm_e)+(e12*f00)*(e12*f00)+(e13*f0n)*(e13*f0n));
    	err_xn0n_c->SetPointEYlow(i,error_x);
    	f_xn0n_eyl[i-1]=error_x;
    	
    	error_x=sqrt((e20*ps_ep)*(e20*ps_ep)+(e21*pm_e)*(e21*pm_e)+(e22*f00)*(e22*f00)+(e23*f0n)*(e23*f0n)+(e24*fnn)*(e24*fnn));
    	//cout << "error_x_high2: " << error_x << endl;
    	err_xnxn_c->SetPointEYhigh(i,error_x);
    	f_xnxn_eyh[i-1]=error_x;
    	error_x=sqrt((e20*ps_en)*(e20*ps_en)+(e21*pm_e)*(e21*pm_e)+(e22*f00)*(e22*f00)+(e23*f0n)*(e23*f0n)+(e24*fnn)*(e24*fnn));
    	err_xnxn_c->SetPointEYlow(i,error_x);
    	f_xnxn_eyl[i-1]=error_x;
    	//cout << "error_x_low2: " << error_x << endl;
    	
    
    	
    }
    
    double fex[2]={0};

    TGraphAsymmErrors* of_0n0n_h = new TGraphAsymmErrors(2,of_0n0n_x,of_0n0n_y,fex,fex,of_0n0n_eyl,of_0n0n_eyh);
    TGraphAsymmErrors* of_xn0n_h = new TGraphAsymmErrors(2,of_xn0n_x,of_xn0n_y,fex,fex,of_xn0n_eyl,of_xn0n_eyh);
    TGraphAsymmErrors* of_xnxn_h = new TGraphAsymmErrors(2,of_xnxn_x,of_xnxn_y,fex,fex,of_xnxn_eyl,of_xnxn_eyh);
    TGraphAsymmErrors* f_0n0n_h = new TGraphAsymmErrors(2,f_0n0n_x,f_0n0n_y,fex,fex,f_0n0n_eyl,f_0n0n_eyh);
    TGraphAsymmErrors* f_xn0n_h = new TGraphAsymmErrors(2,f_xn0n_x,f_xn0n_y,fex,fex,f_xn0n_eyl,f_xn0n_eyh);
    TGraphAsymmErrors* f_xnxn_h = new TGraphAsymmErrors(2,f_xnxn_x,f_xnxn_y,fex,fex,f_xnxn_eyl,f_xnxn_eyh);
    
    cout  << f_0n0n_y[0] << " " << f_0n0n_y[1] << " " << f_0n0n_eyl[1] << " " << f_0n0n_eyh[1] <<  endl;
    cout  << f_xn0n_y[0] << " " << f_xn0n_y[1] << endl;
    cout  << f_xnxn_y[0] << " " << f_xnxn_y[1] << endl;
   
    h_zdc_check_08->Add(h_0n0n08_c,h_xn0n08_c);
    h_zdc_check_08->Add(h_zdc_check_08,h_xnxn08_c);
    
    h_zdc_check_16->Add(h_0n0n16_c,h_xn0n16_c);
    h_zdc_check_16->Add(h_zdc_check_16,h_xnxn16_c);
    
    h_zdc_check_24->Add(h_0n0n24_c,h_xn0n24_c);
    h_zdc_check_24->Add(h_zdc_check_24,h_xnxn24_c);
    
    
    //GetBinContent, SetBinContent, GetBinError, SetBinError
    
    //add MC samples
    h_invm_all->Add(h_ee_invm_high,h_ee_invm_low);
    h_ee_pt_p_all->Add(h_ee_pt_p_high,h_ee_pt_p_low);
    h_ee_pt_n_all->Add(h_ee_pt_n_high,h_ee_pt_n_low);
    h_ee_eta_p_all->Add(h_ee_eta_p_high,h_ee_eta_p_low);
    h_ee_eta_n_all->Add(h_ee_eta_n_high,h_ee_eta_n_low);
    h_ee_phi_p_all->Add(h_ee_phi_p_high,h_ee_phi_p_low);
    h_ee_phi_n_all->Add(h_ee_phi_n_high,h_ee_phi_n_low);
    h_aco_all->Add(h_aco_high,h_aco_low);
    h_pt_ee_all->Add(h_pt_ee_high,h_pt_ee_low);
	h_yy_all->Add(h_yy_high,h_yy_low);
 	
 	h_invm_all->Add(h_invm_all,h_ee_invm_add);
 	h_ee_pt_p_all->Add(h_ee_pt_p_all,h_ee_pt_p_add);
    h_ee_pt_n_all->Add(h_ee_pt_n_all,h_ee_pt_n_add);
    h_ee_eta_p_all->Add(h_ee_eta_p_all,h_ee_eta_p_add);
    h_ee_eta_n_all->Add(h_ee_eta_n_all,h_ee_eta_n_add);
    h_ee_phi_p_all->Add(h_ee_phi_p_all,h_ee_phi_p_add);
    h_ee_phi_n_all->Add(h_ee_phi_n_all,h_ee_phi_n_add);
    h_aco_all->Add(h_aco_all,h_aco_add);
    h_pt_ee_all->Add(h_pt_ee_all,h_pt_ee_add);
    h_yy_all->Add(h_yy_all,h_yy_add);
 
 	//create histograms
 
    TCanvas *c1 = new TCanvas("c1","Electron #eta",700,1500);
    c1->Divide(1,2);
    c1->cd(1)->SetLeftMargin(0.2);
    c1->cd(1)->SetPad(0,0.25,1,1);
    h_ee_eta_p_all->SetTitle("#eta_{e+}");
    h_ee_eta_p_all->SetFillColorAlpha(kBlue-10, 0.35);
    h_ee_eta_p_all->Draw("hist");
    h_ee_eta_p_dat->Draw("same ep");
    h_ee_eta_p_dat->SetMarkerStyle(21);
    h_ee_eta_p_dat->SetLineColor(1);
    h_ee_eta_p_dat->SetMarkerSize(0.5);
    h_ee_eta_p_all->GetXaxis()->SetTitle("#eta");
    h_ee_eta_p_all->GetYaxis()->SetTitle("Counts");
    h_ee_eta_p_all->GetXaxis()->SetRangeUser(-3,3);
    h_ee_eta_p_all->GetYaxis()->SetRangeUser(0,1000);
    TLegend* legend1 = new TLegend(0.230,0.75,0.40,0.85);
    legend1->SetBorderSize(0);
    legend1->SetBorderSize(0);
    legend1->SetFillColor(0);
    legend1->SetTextSize(0.04);    
    legend1->AddEntry(h_ee_eta_p_all,"Starlight","Lf");
    legend1->AddEntry(h_ee_eta_p_dat,"Data","P");
    legend1->Draw();   
    
    c1->cd(2)->SetLeftMargin(0.2);
    c1->cd(2)->SetPad(0,0.0,1,0.25);
    c1->cd(2)->SetBottomMargin(0.25);
    TLine *line_1 = new TLine(-3,1,3,1);
    TH1D* h_ee_eta_p_rat;
    h_ee_eta_p_rat = (TH1D*)h_ee_eta_p_dat->Clone("h_ee_pt_p_rat");
    h_ee_eta_p_rat->GetYaxis()->SetTitle("Data/Starlight");
    h_ee_eta_p_rat->GetXaxis()->SetTitle("#eta");
    h_ee_eta_p_rat->GetXaxis()->SetLabelSize(0.1);
    h_ee_eta_p_rat->GetYaxis()->SetLabelSize(0.1);
    h_ee_eta_p_rat->GetXaxis()->SetTitleSize(0.1);
    h_ee_eta_p_rat->GetYaxis()->SetTitleSize(0.1); //setndivisions
    h_ee_eta_p_rat->GetXaxis()->SetLabelSize(0.1);
    h_ee_eta_p_rat->GetXaxis()->SetTitleOffset(1);
    h_ee_eta_p_rat->GetYaxis()->SetTitleOffset(0.5);
    h_ee_eta_p_rat->GetXaxis()->SetRangeUser(-3,3);
    h_ee_eta_p_rat->SetTitle("");
    h_ee_eta_p_rat->Divide(h_ee_eta_p_all);
    h_ee_eta_p_rat->GetYaxis()->SetRangeUser(0.6,1.4);
    h_ee_eta_p_rat->Draw("ep");
    line_1->Draw("same");
      
    TCanvas *c2 = new TCanvas("c2","Electron #eta",700,1500);
    c2->Divide(1,2);
    c2->cd(1)->SetLeftMargin(0.2);
    c2->cd(1)->SetPad(0,0.25,1,1);
    h_ee_eta_n_all->SetTitle("#eta_{e-}");
    h_ee_eta_n_all->Draw("hist");
    h_ee_eta_n_dat->Draw("same ep");
    h_ee_eta_n_dat->SetLineColor(1);
    h_ee_eta_n_dat->SetMarkerStyle(21);
    h_ee_eta_n_dat->SetMarkerSize(0.5);
    h_ee_eta_n_all->SetFillColorAlpha(kBlue-10, 0.35);
    h_ee_eta_n_all->GetXaxis()->SetRangeUser(-3,3);
    h_ee_eta_n_all->GetYaxis()->SetRangeUser(0,1000);
    h_ee_eta_n_all->GetXaxis()->SetTitle("#eta");
    h_ee_eta_n_all->GetYaxis()->SetTitle("Counts");
    
    
    TLegend* legend2 = new TLegend(0.230,0.75,0.40,0.85);
    legend2->SetBorderSize(0);
    legend2->SetBorderSize(0);
    legend2->SetFillColor(0);
    legend2->SetTextSize(0.04);    
    legend2->AddEntry(h_ee_eta_p_all,"Starlight","Lf");
    legend2->AddEntry(h_ee_eta_p_dat,"Data","P");
    legend2->Draw();
    
    c2->cd(2)->SetLeftMargin(0.2);
    c2->cd(2)->SetPad(0,0.0,1,0.25);
    c2->cd(2)->SetBottomMargin(0.25);
    TLine *line_2 = new TLine(-3,1,3,1);
    TH1D* h_ee_eta_n_rat;
    h_ee_eta_n_rat = (TH1D*)h_ee_eta_n_dat->Clone("h_ee_eta_n_rat");
   	h_ee_eta_n_rat->GetYaxis()->SetTitle("Data/Starlight");
    h_ee_eta_n_rat->GetXaxis()->SetTitle("#eta");
    h_ee_eta_n_rat->GetXaxis()->SetRangeUser(-3,3);
    h_ee_eta_n_rat->Divide(h_ee_eta_n_all);
    h_ee_eta_n_rat->SetTitle("");
    h_ee_eta_n_rat->GetXaxis()->SetLabelSize(0.1);
    h_ee_eta_n_rat->GetYaxis()->SetLabelSize(0.1);
    h_ee_eta_n_rat->GetXaxis()->SetTitleSize(0.1);
    h_ee_eta_n_rat->GetYaxis()->SetTitleSize(0.1); //setndivisions
    h_ee_eta_n_rat->GetXaxis()->SetTitleOffset(1);
    h_ee_eta_n_rat->GetYaxis()->SetTitleOffset(0.5);
    h_ee_eta_n_rat->Draw("ep");
    h_ee_eta_n_rat->GetYaxis()->SetRangeUser(0.6,1.4);
    line_2->Draw("same");
    
    TCanvas *c3 = new TCanvas("c3","Electron pt p",700,1500);
    c3->Divide(1,2);
    c3->cd(1)->SetLeftMargin(0.2);
    c3->cd(1)->SetPad(0,0.25,1,1);
    h_ee_pt_p_all->Draw("hist");
    h_ee_pt_p_dat->Draw("same ep");
    h_ee_pt_p_dat->SetLineColor(1);
    h_ee_pt_p_all->SetTitle("P_{t,e+}");
    h_ee_pt_p_dat->SetMarkerSize(0.5);
    h_ee_pt_p_all->SetFillColorAlpha(kBlue-10, 0.35);
    h_ee_pt_p_dat->SetMarkerStyle(21);
    h_ee_pt_p_all->GetXaxis()->SetTitle("P_{t} [GeV]");
    h_ee_pt_p_all->GetYaxis()->SetTitle("Counts");
    h_ee_pt_p_all->GetYaxis()->SetRangeUser(0,3000);
    TLegend* legend3 = new TLegend(0.65,0.7,0.85,0.85);
    legend3->SetBorderSize(0);
    legend3->SetBorderSize(0);
    legend3->SetFillColor(0);
    legend3->SetTextSize(0.04);    
    legend3->AddEntry(h_ee_eta_p_all,"Starlight","fL");
    legend3->AddEntry(h_ee_eta_p_dat,"Data","P");
    legend3->Draw();
    
    c3->cd(2)->SetLeftMargin(0.2);
    c3->cd(2)->SetPad(0,0.0,1,0.25);
    c3->cd(2)->SetBottomMargin(0.25);
    TLine *line_3 = new TLine(0,1,30,1);
    TH1D* h_ee_pt_p_rat;
    h_ee_pt_p_rat = (TH1D*)h_ee_pt_p_dat->Clone("h_ee_pt_p_rat");
   	h_ee_pt_p_rat->GetYaxis()->SetTitle("Data/Starlight");
    h_ee_pt_p_rat->GetXaxis()->SetTitle("P_{t} [GeV]");
    h_ee_pt_p_rat->Divide(h_ee_pt_p_all);
    h_ee_pt_p_rat->SetTitle("");
    h_ee_pt_p_rat->GetYaxis()->SetRangeUser(0,2);
    h_ee_pt_p_rat->GetXaxis()->SetLabelSize(0.1);
    h_ee_pt_p_rat->GetYaxis()->SetLabelSize(0.1);
    h_ee_pt_p_rat->GetXaxis()->SetTitleSize(0.1);
    h_ee_pt_p_rat->GetYaxis()->SetTitleSize(0.1); //setndivisions
    h_ee_pt_p_rat->GetXaxis()->SetTitleOffset(1);
    h_ee_pt_p_rat->GetYaxis()->SetTitleOffset(0.5);
    h_ee_pt_p_rat->Draw("ep");
    line_3->Draw("same");
    
    TCanvas *c4 = new TCanvas("c4","Electron pt n",700,1500);
    c4->Divide(1,2);
    c4->cd(1)->SetLeftMargin(0.2);
    c4->cd(1)->SetPad(0,0.25,1,1);
    h_ee_pt_n_all->Draw("hist");
    h_ee_pt_n_dat->Draw("same ep");
    h_ee_pt_n_dat->SetLineColor(1);
    h_ee_pt_n_all->SetTitle("P_{t,e-}");
    h_ee_pt_n_all->GetYaxis()->SetRangeUser(0,3000);
    h_ee_pt_n_dat->SetMarkerSize(0.5);
    h_ee_pt_n_all->SetFillColorAlpha(kBlue-10, 0.35);
    h_ee_pt_n_dat->SetMarkerStyle(21);
    h_ee_pt_n_all->GetXaxis()->SetTitle("P_{t} [GeV]");
    h_ee_pt_n_all->GetYaxis()->SetTitle("Counts");
    
    TLegend* legend4 = new TLegend(0.65,0.7,0.85,0.85);
    legend4->SetBorderSize(0);
    legend4->SetBorderSize(0);
    legend4->SetFillColor(0);
    legend4->SetTextSize(0.04);    
    legend4->AddEntry(h_ee_eta_p_all,"Starlight","fL");
    legend4->AddEntry(h_ee_eta_p_dat,"Data","P");
    legend4->Draw();
    
    c4->cd(2)->SetLeftMargin(0.2);
    c4->cd(2)->SetPad(0,0.0,1,0.25);
    c4->cd(2)->SetBottomMargin(0.25);
    TLine *line_4 = new TLine(0,1,30,1);
    TH1D* h_ee_pt_n_rat;
    h_ee_pt_n_rat = (TH1D*)h_ee_pt_n_dat->Clone("h_ee_pt_n_rat");
   	h_ee_pt_n_rat->GetYaxis()->SetTitle("Data/Starlight");
    h_ee_pt_n_rat->GetXaxis()->SetTitle("P_{t} [GeV]");
    h_ee_pt_n_rat->Divide(h_ee_pt_n_all);
    h_ee_pt_n_rat->Draw("ep");
    h_ee_pt_n_rat->SetTitle("");
    h_ee_pt_n_rat->GetYaxis()->SetRangeUser(0,2);
    h_ee_pt_n_rat->GetXaxis()->SetLabelSize(0.1);
    h_ee_pt_n_rat->GetYaxis()->SetLabelSize(0.1);
    h_ee_pt_n_rat->GetXaxis()->SetTitleSize(0.1);
    h_ee_pt_n_rat->GetYaxis()->SetTitleSize(0.1); //setndivisions
    h_ee_pt_n_rat->GetXaxis()->SetTitleOffset(1);
    h_ee_pt_n_rat->GetYaxis()->SetTitleOffset(0.5);
    line_4->Draw("same");
    
    TCanvas *c5 = new TCanvas("c5","Invariant mass",700,1500);
    c5->Divide(1,2);
    c5->cd(1)->SetLeftMargin(0.2);
    c5->cd(1)->SetPad(0,0.25,1,1);
    h_invm_all->Draw("hist");
    h_ee_invm_dat->Draw("same ep");
    h_ee_invm_dat->SetLineColor(1);
    h_ee_invm_dat->SetMarkerStyle(21);
    h_invm_all->SetFillColorAlpha(kBlue-10, 0.35);
    h_invm_all->SetTitle("m_{inv}");
    h_invm_all->GetYaxis()->SetRangeUser(0,3500);
    h_ee_invm_dat->SetMarkerSize(0.5);
    h_invm_all->GetXaxis()->SetTitle("Mass [GeV]");
    h_invm_all->GetYaxis()->SetTitle("Counts");
    
    TLegend* legend5 = new TLegend(0.65,0.7,0.85,0.85);
    legend5->SetBorderSize(0);
    legend5->SetBorderSize(0);
    legend5->SetFillColor(0);
    legend5->SetTextSize(0.04);    
    legend5->AddEntry(h_ee_eta_p_all,"Starlight","fL");
    legend5->AddEntry(h_ee_eta_p_dat,"Data","P");
    legend5->Draw();
    
	c5->cd(2)->SetLeftMargin(0.2);
	c5->cd(2)->SetPad(0,0.0,1,0.25);
	c5->cd(2)->SetBottomMargin(0.25);
    TLine *line_5 = new TLine(0,1,70,1);
    TH1D* h_ee_invm_rat;
    h_ee_invm_rat = (TH1D*)h_ee_invm_dat->Clone("h_ee_invm_rat");
   	h_ee_invm_rat->GetYaxis()->SetTitle("Data/Starlight");
    h_ee_invm_rat->GetXaxis()->SetTitle("Mass [GeV]");
    h_ee_invm_rat->Divide(h_invm_all);
    h_ee_invm_rat->SetTitle("");
    h_ee_invm_rat->GetYaxis()->SetRangeUser(0,2);
    h_ee_invm_rat->GetYaxis()->SetRangeUser(0,2);
    h_ee_invm_rat->GetXaxis()->SetLabelSize(0.1);
    h_ee_invm_rat->GetYaxis()->SetLabelSize(0.1);
    h_ee_invm_rat->GetXaxis()->SetTitleSize(0.1);
    h_ee_invm_rat->GetYaxis()->SetTitleSize(0.1); //setndivisions
    h_ee_invm_rat->GetXaxis()->SetTitleOffset(1);
    h_ee_invm_rat->GetYaxis()->SetTitleOffset(0.5);
    h_ee_invm_rat->Draw("ep");
    line_5->Draw("same");
       
    TCanvas *c6 = new TCanvas("c6","Acoplanarity",700,1500);
    c6->Divide(1,2);
    c6->cd(1)->SetLeftMargin(0.2);
    c6->cd(1)->SetPad(0,0.25,1,1);
    
    h_aco_all->Draw("hist");
    h_aco_dat->Draw("same ep");
    h_aco_dat->SetLineColor(1);
    h_aco_all->SetFillColorAlpha(kBlue-10, 0.35);
    h_aco_all->SetTitle("A");
    h_aco_all->GetYaxis()->SetRangeUser(0,1850);
    h_aco_dat->SetMarkerSize(0.5);
    h_aco_dat->SetMarkerStyle(21);
    h_aco_all->GetXaxis()->SetTitle("Acoplanarity");
    h_aco_all->GetYaxis()->SetTitle("Counts");
    
    TLegend* legend6 = new TLegend(0.40,0.7,0.55,0.85);
    legend6->SetBorderSize(0);
    legend6->SetBorderSize(0);
    legend6->SetFillColor(0);
    legend6->SetTextSize(0.04);    
    legend6->AddEntry(h_ee_eta_p_all,"Starlight","fL");
    legend6->AddEntry(h_ee_eta_p_dat,"Data","P");
    legend6->Draw();  
    
   	c6->cd(2)->SetLeftMargin(0.2);
   	c6->cd(2)->SetPad(0,0.0,1,0.25);
   	c6->cd(2)->SetBottomMargin(0.25);
    TLine *line_6 = new TLine(0,1,0.01,1);
    TH1D* h_aco_rat;
    h_aco_rat = (TH1D*)h_aco_dat->Clone("h_aco_rat");
   	h_aco_rat->GetYaxis()->SetTitle("Data/Starlight");
    h_aco_rat->GetXaxis()->SetTitle("Acoplanarity");
    h_aco_rat->SetTitle("");
    h_aco_rat->GetYaxis()->SetRangeUser(0,2);
    h_aco_rat->GetXaxis()->SetLabelSize(0.1);
    h_aco_rat->GetYaxis()->SetLabelSize(0.1);
    h_aco_rat->GetXaxis()->SetTitleSize(0.1);
    h_aco_rat->GetYaxis()->SetTitleSize(0.1); //setndivisions
    h_aco_rat->GetXaxis()->SetTitleOffset(1);
    h_aco_rat->GetYaxis()->SetTitleOffset(0.5);
    h_aco_rat->GetYaxis()->SetRangeUser(0.6,1.8);
   	h_aco_rat->Divide(h_aco_all);
    h_aco_rat->Draw("ep");
    line_6->Draw("same");
    
    TCanvas *c7 = new TCanvas("c7","Di-electron transverse momentum",700,1500);
    c7->Divide(1,2);
    c7->cd(1)->SetLeftMargin(0.2);
    c7->cd(1)->SetPad(0,0.25,1,1);
    h_pt_ee_all->Draw("hist");
    h_pt_ee_all->GetXaxis()->SetRangeUser(0,2.5);
    h_pt_ee_dat->Draw("same ep");
    h_pt_ee_dat->SetLineColor(1);
    h_pt_ee_dat->SetMarkerStyle(21);
    h_pt_ee_all->SetFillColorAlpha(kBlue-10, 0.35);
    h_pt_ee_dat->SetTitle("p_{t}^{ee}");
    h_pt_ee_dat->GetYaxis()->SetRangeUser(0,1850);
    h_pt_ee_dat->SetMarkerSize(0.5);
    h_pt_ee_all->GetXaxis()->SetTitle("P_{t}^{ee} [GeV]");
    h_pt_ee_all->GetYaxis()->SetTitle("Counts");
    h_pt_ee_all->SetTitle("Di-electron transverse momentum");
    
    TLegend* legend7 = new TLegend(0.65,0.7,0.40,0.85);
    legend7->SetBorderSize(0);
    legend7->SetBorderSize(0);
    legend7->SetFillColor(0);
    legend7->SetTextSize(0.04);    
    legend7->AddEntry(h_ee_eta_p_all,"Starlight","fL");
    legend7->AddEntry(h_ee_eta_p_dat,"Data","P");
    legend7->Draw();  
        
    c7->cd(2)->SetLeftMargin(0.2);
    c7->cd(2)->SetPad(0,0.0,1,0.25);
    c7->cd(2)->SetBottomMargin(0.25);
    TLine *line_7 = new TLine(0,1,2.5,1);
    TH1D* h_pt_ee_rat;
   	h_pt_ee_rat = (TH1D*)h_pt_ee_dat->Clone("h_pt_ee_rat");
   	h_pt_ee_rat->SetTitle("Di-electron transverse momentum");
   	h_pt_ee_rat->GetXaxis()->SetRangeUser(0,2.5);
   	h_pt_ee_rat->GetYaxis()->SetRangeUser(0,2);
   	h_pt_ee_rat->SetTitle("");
   	h_pt_ee_rat->GetYaxis()->SetTitle("Data/Starlight");
    h_pt_ee_rat->GetXaxis()->SetTitle("P_{t}^{ee} [GeV]");
    h_pt_ee_rat->GetXaxis()->SetLabelSize(0.1);
    h_pt_ee_rat->GetYaxis()->SetLabelSize(0.1);
    h_pt_ee_rat->GetXaxis()->SetTitleSize(0.1);
    h_pt_ee_rat->GetYaxis()->SetTitleSize(0.1); //setndivisions
    h_pt_ee_rat->GetXaxis()->SetTitleOffset(1);
    h_pt_ee_rat->GetYaxis()->SetTitleOffset(0.5);
    h_pt_ee_rat->GetYaxis()->SetRangeUser(0.95,1.8);
   	h_pt_ee_rat->Divide(h_pt_ee_all);
    h_pt_ee_rat->Draw("ep");
    line_7->Draw("same");
    
    TCanvas *c8 = new TCanvas("c8","ZDC",1200,800);
    c8->cd(1)->SetLeftMargin(0.2);
    hh_zdc_ee->Draw("LEGO2 0");
    hh_zdc_ee->SetContour(99);
    hh_zdc_ee->GetXaxis()->SetTitle("E_{ZDC_a} / (2.51 TeV)");
    hh_zdc_ee->GetYaxis()->SetTitle("E_{ZDC_c} / (2.51 TeV)");
    hh_zdc_ee->GetXaxis()->SetTitleOffset(2);
    hh_zdc_ee->GetYaxis()->SetTitleOffset(2);
	gPad->SetLogz();
	gPad->SetPhi(225);
    gStyle->SetPalette(kBird);
    hh_zdc_ee->GetXaxis()->SetRangeUser(0,4);
    hh_zdc_ee->GetYaxis()->SetRangeUser(0,4);
    hh_zdc_ee->GetZaxis()->SetRangeUser(0.5,5e4);
    
    TCanvas *c9 = new TCanvas("c9","xn0n,xnxn",800,600);
    
    n_0n0n08->GetXaxis()->SetTitle("m_{ee} [GeV]");
    n_0n0n08->SetTitle("|y_{ee}| < 0.8 ");
    n_0n0n08->GetYaxis()->SetTitle("#it{f_{0n0n}}, #it{f_{Xn0n}}, #it{f_{XnXn}}");
    n_0n0n08->Draw("AP");
    n_xn0n08->Draw("same P");
  	n_xnxn08->Draw("same P");
  	o_0n0n08_h->Draw("same P");
    o_xn0n08_h->Draw("same P");
    o_xnxn08_h->Draw("same P");
    n_0n0n08->SetMarkerColor(2);
    n_xn0n08->SetMarkerColor(2);
    n_xnxn08->SetMarkerColor(2);
    n_xnxn08->SetMarkerStyle(3);
    n_xn0n08->SetMarkerStyle(4);
    n_0n0n08->SetMarkerStyle(5);
    o_xnxn08_h->SetMarkerStyle(3);
    o_xn0n08_h->SetMarkerStyle(4);
    o_0n0n08_h->SetMarkerStyle(5);
  	n_0n0n08->GetYaxis()->SetRangeUser(0,1);
    TLegend* legend9 = new TLegend(0.55,0.45,0.75,0.65);
    legend9->SetBorderSize(0);
    legend9->SetBorderSize(0);
    legend9->SetFillColor(0);
    legend9->SetTextSize(0.02);  
    legend9->SetHeader("corrected", "C");
    legend9->AddEntry(n_0n0n08,"0n0n","P");
    legend9->AddEntry(n_xn0n08,"xn0n","P");
    legend9->AddEntry(n_xnxn08,"xnxn","P");
    legend9->Draw();  
    
    TLegend* legend10 = new TLegend(0.55,0.65,0.75,0.85);
    legend10->SetBorderSize(0);
    legend10->SetBorderSize(0);
    legend10->SetFillColor(0);
    legend10->SetTextSize(0.02);  
    legend10->SetHeader("uncorrected", "C");
    legend10->AddEntry(o_0n0n08_h,"0n0n","P");
    legend10->AddEntry(o_xn0n08_h,"xn0n","P");
    legend10->AddEntry(o_xnxn08_h,"xnxn","P");
    legend10->Draw(); 
    
    TCanvas *c10 = new TCanvas("c10","xn0n,xnxn",800,600);
    n_0n0n16->GetXaxis()->SetTitle("m_{ee} [GeV]");
    n_0n0n16->SetTitle("0.8 < |y_{ee}| < 1.6");
	n_0n0n16->GetYaxis()->SetTitle("#it{f_{0n0n}}, #it{f_{Xn0n}}, #it{f_{XnXn}}");
    n_0n0n16->Draw("AP");
    n_xn0n16->Draw("same P");
    n_xnxn16->Draw("same P");
    o_0n0n16_h->Draw("same P");
    o_xn0n16_h->Draw("same P");
    o_xnxn16_h->Draw("same P");
    n_0n0n16->SetMarkerColor(2);
    n_xn0n16->SetMarkerColor(2);
    n_xnxn16->SetMarkerColor(2);
    n_xnxn16->SetMarkerStyle(3);
    n_xn0n16->SetMarkerStyle(4);
    n_0n0n16->SetMarkerStyle(5);
    o_xnxn16_h->SetMarkerStyle(3);
    o_xn0n16_h->SetMarkerStyle(4);
    o_0n0n16_h->SetMarkerStyle(5);
    n_0n0n16->GetYaxis()->SetRangeUser(0,1);
    
    TLegend* legend11 = new TLegend(0.55,0.45,0.75,0.65);
    legend11->SetBorderSize(0);
    legend11->SetBorderSize(0);
    legend11->SetFillColor(0);
    legend11->SetTextSize(0.02);  
    legend11->SetHeader("corrected", "C");
    legend11->AddEntry(n_0n0n16,"0n0n","P");
    legend11->AddEntry(n_xn0n16,"xn0n","P");
    legend11->AddEntry(n_xnxn16,"xnxn","P");
    legend11->Draw();  
    
    TLegend* legend15 = new TLegend(0.55,0.65,0.75,0.85);
    legend15->SetBorderSize(0);
    legend15->SetBorderSize(0);
    legend15->SetFillColor(0);
    legend15->SetTextSize(0.02);  
    legend15->SetHeader("uncorrected", "C");
    legend15->AddEntry(o_0n0n16_h,"0n0n","P");
    legend15->AddEntry(o_xn0n16_h,"xn0n","P");
    legend15->AddEntry(o_xnxn16_h,"xnxn","P");
    legend15->Draw(); 
    
    TCanvas *c11 = new TCanvas("c11","xn0n,xnxn",800,600);
    n_0n0n24->GetXaxis()->SetTitle("m_{ee} [GeV]");
   	n_0n0n24->SetTitle("1.6 < |y_{ee}| < 2.4");
	n_0n0n24->GetYaxis()->SetTitle("#it{f_{0n0n}}, #it{f_{Xn0n}}, #it{f_{XnXn}}");
    n_0n0n24->Draw("AP");
    n_xn0n24->Draw("same P");
    n_xnxn24->Draw("same P");
    o_0n0n24_h->Draw("same P");
    o_xn0n24_h->Draw("same P");
    o_xnxn24_h->Draw("same P");
    n_0n0n24->SetMarkerColor(2);
    n_xn0n24->SetMarkerColor(2);
    n_xnxn24->SetMarkerColor(2);
   	n_xnxn24->SetMarkerStyle(3);
    n_xn0n24->SetMarkerStyle(4);
    n_0n0n24->SetMarkerStyle(5);
   	o_xnxn24_h->SetMarkerStyle(3);
    o_xn0n24_h->SetMarkerStyle(4);
    o_0n0n24_h->SetMarkerStyle(5);	 
    n_0n0n24->GetYaxis()->SetRangeUser(0,1);
    
    
    
    TLegend* legend16 = new TLegend(0.55,0.45,0.75,0.65);
    legend16->SetBorderSize(0);
    legend16->SetBorderSize(0);
    legend16->SetFillColor(0);
    legend16->SetTextSize(0.02);  
    legend16->SetHeader("corrected", "C");
    legend16->AddEntry(n_0n0n24,"0n0n","P");
    legend16->AddEntry(n_xn0n24,"xn0n","P");
    legend16->AddEntry(n_xnxn24,"xnxn","P");
    legend16->Draw();  
    
    TLegend* legend17 = new TLegend(0.55,0.65,0.75,0.85);
    legend17->SetBorderSize(0);
    legend17->SetBorderSize(0);
    legend17->SetFillColor(0);
    legend17->SetTextSize(0.02);  
    legend17->SetHeader("uncorrected", "C");
    legend17->AddEntry(o_0n0n24_h,"0n0n","P");
    legend17->AddEntry(o_xn0n24_h,"xn0n","P");
    legend17->AddEntry(o_xnxn24_h,"xnxn","P");
    legend17->Draw(); 
    /*
    TCanvas *c9 = new TCanvas("c9","xn0n,xnxn",800,600);
    h_xnxn08->SetMarkerStyle(3);
    //err_0n0n08->Draw("ALP");
    h_xnxn08->Draw("same P");  
    h_xnxn08->GetYaxis()->SetRangeUser(0,1);
    h_xn0n08->Draw("same P");
    h_xn0n08->SetMarkerStyle(4);
    h_0n0n08->Draw("same P");
    h_0n0n08->SetMarkerStyle(5);
    h_0n0n08_c->Draw("same P");
    h_0n0n08_c->SetMarkerStyle(5);
    h_0n0n08_c->SetMarkerColor(2);
    h_xn0n08_c->Draw("same P");
    h_xn0n08_c->SetMarkerStyle(4);
    h_xn0n08_c->SetMarkerColor(2);
    h_xnxn08_c->Draw("same P");
    h_xnxn08_c->SetMarkerStyle(3);
    h_xnxn08_c->SetMarkerColor(2);
	c9->SetLogx();
	h_xnxn08->GetXaxis()->SetTitle("m_{ee} [GeV]");
    h_xnxn08->GetYaxis()->SetTitle("f_{xn0n}, f_{xnxn}");
	
	TLegend* legend9 = new TLegend(0.25,0.70,0.45,0.85);
    legend9->SetBorderSize(0);
    legend9->SetBorderSize(0);
    legend9->SetFillColor(0);
    legend9->SetTextSize(0.04);  
    legend9->SetHeader("|y_{ee}| < 0.8", "C");
    legend9->AddEntry(h_0n0n08,"0n0n","P");
    legend9->AddEntry(h_xn0n08,"xn0n","P");
    legend9->AddEntry(h_xnxn08,"xnxn","P");
    legend9->Draw();  
    
    TCanvas *c10 = new TCanvas("c10","xn0n,xnxn",800,600);
    h_xnxn16->SetMarkerStyle(3);
    h_xnxn16->Draw("eP");
    h_xnxn16->GetYaxis()->SetRangeUser(0,1);
    h_xn0n16->Draw("same eP");
    h_xn0n16->SetMarkerStyle(4);
    h_0n0n16->Draw("same eP");
    h_0n0n16->SetMarkerStyle(5);
    h_0n0n16_c->Draw("same eP");
    h_0n0n16_c->SetMarkerStyle(5);
    h_0n0n16_c->SetMarkerColor(2);
    h_xn0n16_c->Draw("same eP");
    h_xn0n16_c->SetMarkerStyle(4);
    h_xn0n16_c->SetMarkerColor(2);
    h_xnxn16_c->Draw("same eP");
    h_xnxn16_c->SetMarkerStyle(3);
    h_xnxn16_c->SetMarkerColor(2);
	c10->SetLogx();
	h_xnxn16->GetXaxis()->SetTitle("m_{ee} [GeV]");
    h_xnxn16->GetYaxis()->SetTitle("f_{xn0n}, f_{xnxn}");
	
	TLegend* legend10 = new TLegend(0.25,0.70,0.45,0.85);
    legend10->SetBorderSize(0);
    legend10->SetBorderSize(0);
    legend10->SetFillColor(0);
    legend10->SetTextSize(0.04);  
    legend10->SetHeader("0.8 < |y_{ee}| < 1.6", "C");
    legend10->AddEntry(h_0n0n16,"0n0n","P");
    legend10->AddEntry(h_xn0n16,"xn0n","P");
    legend10->AddEntry(h_xnxn16,"xnxn","P");
    legend10->Draw();  
    
    TCanvas *c11 = new TCanvas("c11","xn0n,xnxn",800,600);
    h_xnxn24->SetMarkerStyle(3);
    h_xnxn24->Draw("eP");
    h_xnxn24->GetYaxis()->SetRangeUser(0,1.0);
    h_xn0n24->Draw("same eP");
    h_xn0n24->SetMarkerStyle(4);
    h_0n0n24->Draw("same eP");
    h_0n0n24->SetMarkerStyle(5);
    h_0n0n24_c->Draw("same eP");
    h_0n0n24_c->SetMarkerStyle(5);
    h_0n0n24_c->SetMarkerColor(2);
    h_xn0n24_c->Draw("same eP");
    h_xn0n24_c->SetMarkerStyle(4);
    h_xn0n24_c->SetMarkerColor(2);
    h_xnxn24_c->Draw("same eP");
    h_xnxn24_c->SetMarkerStyle(3);
    h_xnxn24_c->SetMarkerColor(2);

	c11->SetLogx();
	h_xnxn24->GetXaxis()->SetTitle("m_{ee} [GeV]");
    h_xnxn24->GetYaxis()->SetTitle("f_{xn0n}, f_{xnxn}");
	
	TLegend* legend11 = new TLegend(0.25,0.70,0.45,0.85);
    legend11->SetBorderSize(0);
    legend11->SetBorderSize(0);
    legend11->SetFillColor(0);
    legend11->SetTextSize(0.04);  
    legend11->SetHeader("1.6 < |y_{ee}| < 2.4", "C");
    legend11->AddEntry(h_0n0n24,"0n0n","P");
    legend11->AddEntry(h_xn0n24,"xn0n","P");
    legend11->AddEntry(h_xnxn24,"xnxn","P");
    legend11->Draw(); 
   */
   TCanvas *c12 = new TCanvas("c12","Xn0n,XnXn",800,600);
   	of_0n0n_h->Draw("AP");
   	of_xn0n_h->Draw("same P");
   	of_xnxn_h->Draw("same P");
   	f_0n0n_h->Draw("same P");
   	f_xn0n_h->Draw("same P");
   	f_xnxn_h->Draw("same P");
   	//h_f0n0n->GetXaxis()->SetBinLabel(1,"p_{T}^{e}>2.5 GeV, M_{inv}>6 GeV");  
   	//h_f0n0n->GetXaxis()->SetBinLabel(2,"p_{T}^{e}>4 GeV, M_{inv}>10 GeV"); 
   	of_0n0n_h->GetYaxis()->SetTitle("#it{f_{0n0n}}, #it{f_{Xn0n}}, #it{f_{XnXn}}");
   	of_0n0n_h->SetTitle("");
	f_0n0n_h->SetMarkerColor(2);
   	f_xn0n_h->SetMarkerColor(2);
   	f_xnxn_h->SetMarkerColor(2);
   	of_0n0n_h->GetYaxis()->SetRangeUser(0,1);
   	of_0n0n_h->GetXaxis()->SetLimits(0.0,2.0);
   	of_0n0n_h->SetMarkerStyle(4);
   	of_xn0n_h->SetMarkerStyle(3);
   	of_xnxn_h->SetMarkerStyle(2);
   	
   	f_0n0n_h->SetMarkerStyle(4);
   	f_xn0n_h->SetMarkerStyle(3);
   	f_xnxn_h->SetMarkerStyle(2);
	
	TLegend* legend12 = new TLegend(0.15,0.75,0.35,0.85);
	legend12->SetHeader("0n0n");
	legend12->SetBorderSize(0);
	legend12->SetBorderSize(0);
	legend12->SetFillStyle(0);
	legend12->SetTextSize(0.03);
	legend12->AddEntry(of_0n0n_h,"Uncorrected","p");
	legend12->AddEntry(f_0n0n_h,"Pileup corrected","p");
	legend12->Draw();
	
	TLegend* legend13 = new TLegend(0.40,0.75,0.60,0.85);
	legend13->SetHeader("Xn0n");
	legend13->SetBorderSize(0);
	legend13->SetBorderSize(0);
	legend13->SetFillStyle(0);
	legend13->SetTextSize(0.03);
	legend13->AddEntry(of_xn0n_h,"Uncorrected","p");
	legend13->AddEntry(f_xn0n_h,"Pileup corrected","p");
	legend13->Draw();
	
	TLegend* legend14 = new TLegend(0.65,0.75,0.9,0.85);
	legend14->SetHeader("XnXn");
	legend14->SetBorderSize(0);
	legend14->SetBorderSize(0);
	legend14->SetFillStyle(0);
	legend14->SetTextSize(0.03);
	legend14->AddEntry(of_xnxn_h,"Uncorrected","p");
	legend14->AddEntry(f_xn0n_h,"Pileup corrected","p");
	legend14->Draw();
	
	TCanvas *c20 = new TCanvas("c20","#phi_{+}",700,1500);
    c20->Divide(1,2);
    c20->cd(1)->SetLeftMargin(0.2);
    c20->cd(1)->SetPad(0,0.25,1,1);
    h_ee_phi_p_all->Draw("hist");
    h_ee_phi_p_dat->Draw("same ep");
    h_ee_phi_p_dat->SetLineColor(1);
    h_ee_phi_p_all->SetTitle("#phi_{e+}");
    h_ee_phi_p_dat->SetMarkerSize(0.5);
    h_ee_phi_p_all->SetFillColorAlpha(kBlue-10, 0.35);
    h_ee_phi_p_dat->SetMarkerStyle(21);
    h_ee_phi_p_all->GetXaxis()->SetTitle("#phi_{+} [rad]");
    h_ee_phi_p_all->GetYaxis()->SetTitle("Counts");
    h_ee_phi_p_all->GetYaxis()->SetRangeUser(0,800);
    h_ee_phi_p_all->GetXaxis()->SetRangeUser(-3.5,3.5);
    TLegend* legend20 = new TLegend(0.65,0.7,0.85,0.85);
    legend20->SetBorderSize(0);
    legend20->SetBorderSize(0);
    legend20->SetFillColor(0);
    legend20->SetTextSize(0.04);    
    legend20->AddEntry(h_ee_eta_p_all,"Starlight","fL");
    legend20->AddEntry(h_ee_eta_p_dat,"Data","P");
    legend20->Draw();
    
    c20->cd(2)->SetLeftMargin(0.2);
    c20->cd(2)->SetPad(0,0.0,1,0.25);
    c20->cd(2)->SetBottomMargin(0.25);
    TLine *line_20 = new TLine(-3.5,1,3.5,1);
    TH1D* h_ee_phi_p_rat;
    h_ee_phi_p_rat = (TH1D*)h_ee_phi_p_dat->Clone("h_ee_phi_p_rat");
   	h_ee_phi_p_rat->GetYaxis()->SetTitle("Data/Starlight");
    h_ee_phi_p_rat->GetXaxis()->SetTitle("#phi_{+} [rad]");
    h_ee_phi_p_rat->Divide(h_ee_phi_p_all);
    h_ee_phi_p_rat->GetXaxis()->SetRangeUser(-3.5,3.5);
    h_ee_phi_p_rat->GetYaxis()->SetRangeUser(0.8,1.5);
    h_ee_phi_p_rat->SetTitle("");
    h_ee_phi_p_rat->GetXaxis()->SetLabelSize(0.1);
    h_ee_phi_p_rat->GetYaxis()->SetLabelSize(0.1);
    h_ee_phi_p_rat->GetXaxis()->SetTitleSize(0.1);
    h_ee_phi_p_rat->GetYaxis()->SetTitleSize(0.1); //setndivisions
    h_ee_phi_p_rat->GetXaxis()->SetTitleOffset(1);
    h_ee_phi_p_rat->GetYaxis()->SetTitleOffset(0.5);
    h_ee_phi_p_rat->Draw("ep");
    line_20->Draw("same");
    
    TCanvas *c21 = new TCanvas("c21","#phi_{-}",700,1500);
    c21->Divide(1,2);
    c21->cd(1)->SetLeftMargin(0.2);
    c21->cd(1)->SetPad(0,0.25,1,1);
    h_ee_phi_n_all->Draw("hist");
    h_ee_phi_n_dat->Draw("same ep");
    h_ee_phi_n_dat->SetLineColor(1);
    h_ee_phi_n_all->SetTitle("#phi_{e-}");
    h_ee_phi_n_dat->SetMarkerSize(0.5);
    h_ee_phi_n_all->SetFillColorAlpha(kBlue-10, 0.35);
    h_ee_phi_n_dat->SetMarkerStyle(21);
    h_ee_phi_n_all->GetXaxis()->SetTitle("#phi_{-} [rad]");
    h_ee_phi_n_all->GetYaxis()->SetTitle("Counts");
    h_ee_phi_n_all->GetYaxis()->SetRangeUser(0,800);
    h_ee_phi_n_all->GetXaxis()->SetRangeUser(-3.5,3.5);
    TLegend* legend21 = new TLegend(0.65,0.7,0.85,0.85);
    legend21->SetBorderSize(0);
    legend21->SetBorderSize(0);
    legend21->SetFillColor(0);
    legend21->SetTextSize(0.04);    
    legend21->AddEntry(h_ee_eta_p_all,"Starlight","fL");
    legend21->AddEntry(h_ee_eta_p_dat,"Data","P");
    legend21->Draw();
    
    c21->cd(2)->SetLeftMargin(0.2);
    c21->cd(2)->SetPad(0,0.0,1,0.25);
    c21->cd(2)->SetBottomMargin(0.25);
    TLine *line_21 = new TLine(-3.5,1,3.5,1);
    TH1D* h_ee_phi_n_rat;
    h_ee_phi_n_rat = (TH1D*)h_ee_phi_n_dat->Clone("h_ee_phi_n_rat");
   	h_ee_phi_n_rat->GetYaxis()->SetTitle("Data/Starlight");
    h_ee_phi_n_rat->GetXaxis()->SetTitle("#phi_{-} [rad]");
    h_ee_phi_n_rat->Divide(h_ee_phi_p_all);
    h_ee_phi_n_rat->GetXaxis()->SetRangeUser(-3.5,3.5);
    h_ee_phi_n_rat->GetYaxis()->SetRangeUser(0.8,1.5);
    h_ee_phi_n_rat->SetTitle("");
    h_ee_phi_n_rat->GetXaxis()->SetLabelSize(0.1);
    h_ee_phi_n_rat->GetYaxis()->SetLabelSize(0.1);
    h_ee_phi_n_rat->GetXaxis()->SetTitleSize(0.1);
    h_ee_phi_n_rat->GetYaxis()->SetTitleSize(0.1); //setndivisions
    h_ee_phi_n_rat->GetXaxis()->SetTitleOffset(1);
    h_ee_phi_n_rat->GetYaxis()->SetTitleOffset(0.5);
    h_ee_phi_n_rat->Draw("ep");
    line_21->Draw("same");
	
	
	TCanvas *c22 = new TCanvas("c22","p_{t1} vs p_{t2}",700,1500);
	hh_pt_data->Draw("scat");
	hh_pt_data->SetTitle("p_{t1} vs p_{t2} data");
	hh_pt_data->GetXaxis()->SetTitle("p_{t1} [GeV]");
	hh_pt_data->GetYaxis()->SetTitle("p_{t2} [GeV]");
	
	TCanvas *c23 = new TCanvas("c23","p_{t1} vs p_{t2}",700,1500);
	hh_pt_low->Draw("scat");
	hh_pt_high->Draw("scat same");
	hh_pt_mid->Draw("scat same");
	hh_pt_low->SetTitle("p_{t1} vs p_{t2} MC");
	hh_pt_low->GetXaxis()->SetTitle("p_{t1} [GeV]");
	hh_pt_low->GetYaxis()->SetTitle("p_{t2} [GeV]");
	
	
   	//TH2D* hh_pt_low = (TH2D*)mc_file_low->Get("hh_pt1_pt2"); 
   	//TH2D* hh_pt_high = (TH2D*)mc_file_high->Get("hh_pt1_pt2");
   	//TH2D* hh_pt_mid = (TH2D*)mc_file_add->Get("hh_pt1_pt2");
	
   /*
   	TCanvas *c12 = new TCanvas("c12","Xn0n,XnXn",800,600);
   	h_f0n0n->Draw("eP");
   	h_fxn0n->Draw("same eP");
   	h_fxnxn->Draw("same eP");
   	h_f0n0n_c->Draw("same eP");
   	h_fxn0n_c->Draw("same eP");
   	h_fxnxn_c->Draw("same eP");
   	h_f0n0n->GetXaxis()->SetBinLabel(1,"p_{T}^{e}>3 GeV, M_{inv}>6 GeV");  
   	h_f0n0n->GetXaxis()->SetBinLabel(2,"p_{T}^{e}>4 GeV, M_{inv}>10 GeV"); 
   	h_f0n0n->GetYaxis()->SetTitle("#it{f_{0n0n}}, #it{f_{Xn0n}}, #it{f_{XnXn}}");
	h_f0n0n_c->SetMarkerColor(2);
   	h_fxn0n_c->SetMarkerColor(2);
   	h_fxnxn_c->SetMarkerColor(2);
   	h_f0n0n->GetYaxis()->SetRangeUser(0,1);
   	h_f0n0n->SetMarkerStyle(4);
   	h_fxn0n->SetMarkerStyle(3);
   	h_fxnxn->SetMarkerStyle(2);
   	
   	h_f0n0n_c->SetMarkerStyle(4);
   	h_fxn0n_c->SetMarkerStyle(3);
   	h_fxnxn_c->SetMarkerStyle(2);
	
	TLegend* legend12 = new TLegend(0.15,0.75,0.35,0.85);
	legend12->SetHeader("0n0n");
	legend12->SetBorderSize(0);
	legend12->SetBorderSize(0);
	legend12->SetFillStyle(0);
	legend12->SetTextSize(0.03);
	legend12->AddEntry(h_f0n0n,"Uncorrected","p");
	legend12->AddEntry(h_f0n0n_c,"Pileup corrected","p");
	legend12->Draw();
	
	TLegend* legend13 = new TLegend(0.40,0.75,0.60,0.85);
	legend13->SetHeader("Xn0n");
	legend13->SetBorderSize(0);
	legend13->SetBorderSize(0);
	legend13->SetFillStyle(0);
	legend13->SetTextSize(0.03);
	legend13->AddEntry(h_fxn0n,"Uncorrected","p");
	legend13->AddEntry(h_fxn0n_c,"Pileup corrected","p");
	legend13->Draw();
	
	TLegend* legend14 = new TLegend(0.65,0.75,0.9,0.85);
	legend14->SetHeader("XnXn");
	legend14->SetBorderSize(0);
	legend14->SetBorderSize(0);
	legend14->SetFillStyle(0);
	legend14->SetTextSize(0.03);
	legend14->AddEntry(h_fxnxn,"Uncorrected","p");
	legend14->AddEntry(h_fxnxn_c,"Pileup corrected","p");
	legend14->Draw();
	*/
    /*
    TCanvas *c13 = new TCanvas("c13","yy_ee",800,600);
    c13->Divide(1,2);
    c13->cd(1)->SetLeftMargin(0.2);
    h_yy_all->Draw("hist");
    h_yy_dat->Draw("same p");
    h_yy_dat->SetMarkerStyle(20);
    h_yy_all->GetXaxis()->SetTitle("yy_{ee}");
    h_yy_all->GetYaxis()->SetTitle("Counts");
    
    TLegend* legend15 = new TLegend(0.25,0.4,0.40,0.65);
    legend15->SetBorderSize(0);
    legend15->SetBorderSize(0);
    legend15->SetFillColor(0);
    legend15->SetTextSize(0.04);    
    legend15->AddEntry(h_ee_eta_p_all,"Starlight","L");
    legend15->AddEntry(h_ee_eta_p_dat,"Data","P");
    legend15->Draw();  
    
    c13->cd(2)->SetLeftMargin(0.2);
    TLine *line_13 = new TLine(-5,1,5,1);
    TH1D* h_yy_rat;
    h_yy_rat = (TH1D*)h_yy_dat->Clone("h_ee_invm_rat");
   	h_yy_rat->GetYaxis()->SetTitle("Data/Starlight");
    h_yy_rat->GetXaxis()->SetTitle("yy_{ee}");
    h_yy_rat->Divide(h_yy_all);
    h_yy_rat->GetYaxis()->SetRangeUser(0,2);
    h_yy_rat->Draw("ep");
    line_13->Draw("same");
    */
    
    //calculate integrals
    
    double sum_high = h_ee_pt_p_high->Integral();
    double sum_low = h_ee_pt_p_low->Integral();
    double sum_data = h_ee_eta_p_dat->Integral();
    
    cout << sum_high << " " << sum_low << " "  << sum_high+sum_low << " " << sum_data << endl;
    
  	//save histograms as images and .root file

    c1->SaveAs("plot_eta_p.png"); 
    c2->SaveAs("plot_eta_n.png");
    c3->SaveAs("plot_pt_p.png");
    c4->SaveAs("plot_pt_n.png");
    c5->SaveAs("plot_invm.png");
    c6->SaveAs("plot_acoll.png");
    c7->SaveAs("plot_pt_al.png");
    c8->SaveAs("plot_zdc.png");
    c9->SaveAs("plot_08_post.png");
    c10->SaveAs("plot_16_post.png");
    c11->SaveAs("plot_24_post.png");
    c12->SaveAs("F_post.png");
    c20->SaveAs("plot_phi_p.png");
    c21->SaveAs("plot_phi_n.png");
    c22->SaveAs("ptpt_dat.png");
    c23->SaveAs("ptpt_mc.png");
    
    TFile *outHistFile =TFile::Open("output_all.root","RECREATE");
    outHistFile->cd();
    h_invm_all->Write();
    h_ee_pt_p_all->Write();
    h_ee_pt_n_all->Write();
    h_ee_eta_p_all->Write();
    h_ee_eta_n_all->Write();
    h_ee_phi_p_all->Write();
    h_ee_phi_n_all->Write();
    h_aco_all->Write();
    outHistFile->Close();
    
    
}
    
