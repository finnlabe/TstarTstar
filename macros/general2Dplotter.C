// Store 2D histograms filled in UHH2 model as pdf
// author: A.Karavdina
// date: 23.09.2019
// Run it with following command:
// root -l -b -q general2Dplotter.C

void general2Dplotter(TString filename="MC.TstarTstar_Combined.root", TString subpath="DNN_Hists", TString histname="2D_DNN_1"){
  // gStyle->SetOptStat(0);
  // gStyle->SetTitleSize(0.045,"x");  
  // gStyle->SetTitleSize(0.045,"y");
  // gStyle->SetTitleYOffset(0.9);

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

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.17);

  gROOT->ForceStyle();
  Double_t w = 600;
  Double_t h = 600;
  

 //Files after selection
 TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/MCStudy/hadded/";
 TString fileprefix = "uhh2.AnalysisModuleRunner.";

 TFile *input = TFile::Open(path+fileprefix+filename);
 TH2D *hist = (TH2D*)input->Get(subpath+"/"+histname);

 TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
 // c1_hist->SetLogz(1);
 //hist->GetZaxis()->SetRangeUser(0,1);
 hist->Draw("colz");
 c1_hist->SaveAs(filename+"_"+subpath+"_"+histname+".pdf");

}
