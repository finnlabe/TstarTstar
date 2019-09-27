// Efficiency histograms for trigger related studies
// author: A.Karavdina
// date: 24.09.2019
// Run it with following command:
// root -l -b -q TriggerEffPlots.C

void TriggerEffPlots(TString filename="uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgammaTgamma_M-700_Run2016v3.root", TString label="TstarTstarToTgammaTgamma_M-700_Run2016v3", TString subpath="SemiLepTTBarMatchGENRECO_triggerSingleLeptonMu",TString histname="Pt_mu", TString subdenomname="SemiLepTTBarMatchGENRECO_mu"){
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
  //We expect histograms filled with and without trigger selection stored in the same file
  //  TString path = "/nfs/dust/cms/user/karavdia/TstarTstar/102X_v1/Preselection/RunII_2016_MuonHihjPtId_cutBasedPhotonIDlooseFall17/semileptonic/";
  TString path = "/nfs/dust/cms/user/karavdia/TstarTstar/102X_v1/Preselection/RunII_2016_MuonHihjPtId_cutBasedPhotonIDlooseFall17_nonIsoHLT/semileptonic/";
  TFile *input = TFile::Open(path+filename);
  TH1D *hist_trigger = (TH1D*)input->Get(subpath+"/"+histname);//histogram after trigger selection
  TH1D *hist_denom = (TH1D*)input->Get(subdenomname+"/"+histname);//histogram before trigger selection
  if(!hist_trigger || !hist_denom) cout<<"Hists are empty"<<endl;;
  if(!hist_trigger || !hist_denom) return;
  if(hist_denom->GetEntries()>0 && hist_trigger->GetEntries()>0) hist_trigger->Divide(hist_denom);
  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  hist_trigger->GetXaxis()->SetTitle(hist_trigger->GetTitle());
  hist_trigger->GetYaxis()->SetTitle("Efficiency");
  hist_trigger->GetYaxis()->SetRangeUser(0,1.5);
  hist_trigger->SetTitle("");
  hist_trigger->SetMarkerStyle(20);
  hist_trigger->SetMarkerColor(1);
  hist_trigger->Draw();
  c1_hist->SaveAs("TrgEff_"+label+"_"+subpath+"_"+histname+".pdf");
}
