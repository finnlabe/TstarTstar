// Efficiency histograms for trigger related studies
// author: F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q cutflowPlots.C

void cutflowPlot(TString filename="uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgammaTgluon_M-1500_Run2016v3.root", TString label="TstarTstarToTgammaTgluon_M-1500_Run2016v3"){
    
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.06,"x");  
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"x");  
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleYOffset(1.20);
  gStyle->SetTitleXOffset(1.0);
  gStyle->SetPalette(55);

  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.17);

  gROOT->ForceStyle();
  Double_t w = 600;
  Double_t h = 600;

  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Preselection/RunII_2016_MuonHihjPtId_mvaPhoIDwp90Fall17_nonIsoandIsoHLT_addTTBarRECO/semileptonic/";
  TString histname = "N_jets";

  int STEPCOUNT = 6;
  TString steps[6] = {"AfterCommon", "AfterNjets", "AfterLepSel", "AfterNpho", "After2D", "After_TstarTstar_Reco"};

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  c1_hist->SetLogy();
  
  TFile *input = TFile::Open(path+filename);
  TH1D *hist_initial = (TH1D*)input->Get(steps[0]+"/"+histname);
  double initial = hist_initial->GetEntries();
  if(initial == 0){initial = 99999;}
  
  TH1* cutflow = new TH1D("cutflow", "cutflow", STEPCOUNT, 0, STEPCOUNT);

  for(int i = 0; i < STEPCOUNT; i++){
    TFile *input = TFile::Open(path+filename);
    TH1D *hist = (TH1D*)input->Get(steps[i]+"/"+histname);
    int entries = hist->GetEntries();
    std::cout << "Entries: " << entries << endl;
    
    cutflow->GetXaxis()->SetBinLabel(i+1,steps[i]);
    cutflow->SetBinContent(i+1, entries/initial);
  }
  
  cutflow->GetYaxis()->SetTitle("Efficiency");
  cutflow->Draw("hist");
  c1_hist->SaveAs("Cutflow_"+label+".pdf");
}
