// Efficiency histograms for trigger related studies
// author: A.Karavdina
// date: 24.09.2019
// Run it with following command:
// root -l -b -q TriggerEffPlots.C

void TriggerEffPlots(TString filename="uhh2.AnalysisModuleRunner.MC.TstarTstar_M-700.root", TString label="TstarTstar_M-700", TString subpath="triggerSingleLeptonEle",TString histname="pt_ele", TString subdenomname="After2D_ele"){
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
  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Selection/hadded/";
  TFile *input = TFile::Open(path+filename);
  TH1D *hist_trigger = (TH1D*)input->Get(subpath+"/"+histname);//histogram after trigger selection
  TH1D *hist_denom = (TH1D*)input->Get(subdenomname+"/"+histname);//histogram before trigger selection

  if(!hist_trigger || !hist_denom) cout<<"Hists are empty"<<endl;;
  if(!hist_trigger || !hist_denom) return;
  //TEfficiency eff;
  TGraphAsymmErrors eff = TGraphAsymmErrors();
  //if(hist_denom->GetEntries()>0 && hist_trigger->GetEntries()>0) eff=TEfficiency(*hist_trigger,*hist_denom);
  if(hist_denom->GetEntries()>0 && hist_trigger->GetEntries()>0) eff.Divide(hist_trigger, hist_denom, "cl=0.68 b(1,1) mode");
  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  eff.GetXaxis()->SetTitle(hist_trigger->GetTitle());
  eff.GetYaxis()->SetTitle("Efficiency");
  //eff.GetYaxis()->SetRangeUser(0,1.5);
  eff.SetTitle("");
  eff.SetMarkerStyle(20);
  eff.SetMarkerColor(1);
  eff.Draw("AP");
  c1_hist->SaveAs("TrgEff_"+label+"_"+subpath+"_"+histname+".pdf");
}
