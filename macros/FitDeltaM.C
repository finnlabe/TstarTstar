// Just a basic macro to create simple plots from given hists.
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q FitDeltaM.C

void FitDeltaM(TString channel = "tgtg", TString masspoint="700", TString subpath="RecoPlots_GEN", TString histname="deltaM_Tstar"){
  
  Double_t w = 800;
  Double_t h = 600;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0001);

  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/MCStudy/tgtg/hadded/";
  TString filename = "uhh2.AnalysisModuleRunner.MC.TstarTstar_M-"+masspoint+".root";
  
  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);

  TFile *input = TFile::Open(path+filename);
  TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname); //histogram
  if(!hist) cout<<"Hist is empty"<<endl;
  hist->Rebin(8);
  TF1 *fit = new TF1("fit", "gaus", -2000, 2000);
  hist->Fit(fit, "", "SAME", -500, 500);
  hist->GetXaxis()->SetTitle("#DeltaM_{T*T*} [GeV]");
  hist->GetYaxis()->SetTitle("Events / GeV");
  hist->SetTitle("");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(1);
  hist->SetLineColor(1);
  hist->Draw("hist");
  fit->Draw("same");

  c1_hist->SaveAs(channel+"_M-"+masspoint+"_"+subpath+"_"+histname+".pdf");
}
