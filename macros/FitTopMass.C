// Just a basic macro to create simple plots from given hists.
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q FitDeltaM.C

void FitTopMass(TString channel = "tgtg", TString subpath="RecoPlots_GEN_onlyttbar_ttag", TString histname="Mtophad"){
  
  Double_t w = 800;
  Double_t h = 600;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0001);
  gStyle->SetHistFillStyle(1);


  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/MCStudy/HOTVR/hadded/";
  TString filename = "FullSignal.root";
  
  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);

  TFile *input = TFile::Open(path+filename);
  TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname); //histogram
  if(!hist) cout<<"Hist is empty"<<endl;
  TF1 *fit = new TF1("fit", "gaus", 0, 500);
  hist->Fit(fit, "", "SAME", 125, 220);
  hist->GetXaxis()->SetTitle("M_{t} [GeV]");
  hist->GetYaxis()->SetTitle("Events / GeV");
  hist->SetTitle("");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(1);
  hist->SetFillColor(kGray);
  hist->SetLineColor(1);
  hist->Draw("hist");
  fit->Draw("same");

  c1_hist->SaveAs("Fit_"+subpath+"_"+histname+".pdf");
}
