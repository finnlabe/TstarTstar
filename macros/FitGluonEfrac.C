// Just a basic macro to create simple plots from given hists.
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q FitDeltaM.C

void FitGluonEfrac(){
  
  Double_t w = 800;
  Double_t h = 600;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0001);
  gStyle->SetHistFillStyle(1);


  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/GenInfo/tgtg/forFits/";
  TString filename = "allSignals.root";
  TString subpath = "GenInfo";
  
  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);

  TFile *input = TFile::Open(path+filename);
  TH1D *hist = (TH1D*)input->Get(subpath+"/"+"gluon_Efrac_AK4"); //histogram
  if(!hist) cout<<"Hist is empty"<<endl;
  TF1 *fit = new TF1("fit", "gaus", 0.5, 1.5);
  //hist->Fit(fit, "", "SAME", 0, 2);
  hist->GetXaxis()->SetTitle("E_{jet}/E_{gluon}");
  hist->GetYaxis()->SetTitle("Events / GeV");
  hist->SetTitle("");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(1);
  hist->SetFillColor(kGray);
  hist->SetLineColor(1);
  hist->Draw("hist");
  //fit->Draw("same");
  c1_hist->SaveAs("Gluon_Efrac_AK4.pdf");

  TH1D *hist2 = (TH1D*)input->Get(subpath+"/"+"gluon_Efrac_AK8"); //histogram
  if(!hist2) cout<<"Hist is empty"<<endl;
  TF1 *fit2 = new TF1("fit", "gaus", 0.5, 1.5);
  //hist2->Fit(fit2, "", "SAME", 0, 2);
  hist2->GetXaxis()->SetTitle("E_{jet}/E_{gluon}");
  hist2->GetYaxis()->SetTitle("Events / GeV");
  hist2->SetTitle("");
  hist2->SetMarkerStyle(20);
  hist2->SetMarkerColor(1);
  hist2->SetFillColor(kGray);
  hist2->SetLineColor(1);
  hist2->Draw("hist");
  //fit2->Draw("same");
  c1_hist->SaveAs("Gluon_Efrac_AK8.pdf");
}
