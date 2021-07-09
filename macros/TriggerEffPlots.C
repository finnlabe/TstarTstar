// Efficiency histograms for trigger related studies
// author: A.Karavdina, changes F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q TriggerEffPlots.C

void TriggerEffPlots(TString filename="uhh2.AnalysisModuleRunner.MC.TstarTstar.root", TString label="TstarTstar"){

  bool bPlotRatio = false;

  Double_t w = 400;
  Double_t h = 400;

  TCanvas *m_can = new TCanvas("canvas", "canvas", w, h);

  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.16);  gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.05);

  // ###################### end cosmetics #######################

  TString channel="ele";
  TString subpath_pre="TriggerXcheck";
  TString subpath_post="AfterST";

  //Files after selection
  //We expect histograms filled with and without trigger selection stored in the same file
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/hadded/";
  TFile *input = TFile::Open(path+filename);
  TH1D *hist_trigger = (TH1D*)input->Get(subpath_post+"_"+channel+"/pt_"+channel);//histogram after trigger selection
  TH1D *hist_denom = (TH1D*)input->Get(subpath_pre+"_"+channel+"/pt_"+channel);//histogram before trigger selection

  hist_trigger->Rebin(5);
  hist_denom->Rebin(5);

  hist_trigger->Scale(0.999);

  if(!hist_trigger || !hist_denom) cout<<"Hists are empty"<<endl;;
  if(!hist_trigger || !hist_denom) return;
  //TEfficiency eff;
  TGraphAsymmErrors eff = TGraphAsymmErrors();
  //if(hist_denom->GetEntries()>0 && hist_trigger->GetEntries()>0) eff=TEfficiency(*hist_trigger,*hist_denom);
  if(hist_denom->GetEntries()>0 && hist_trigger->GetEntries()>0) eff.Divide(hist_trigger, hist_denom, "cl=0.68 b(1,1) mode");

  //eff.Scale(1/0.999);
  for(int i=0;i<eff.GetN();i++) eff.GetY()[i] /= 0.999;

  // cosmetics
  eff.GetXaxis()->SetTitle(hist_trigger->GetTitle());
  eff.GetYaxis()->SetTitle("Efficiency");
  eff.SetTitle("");

  eff.GetXaxis()->SetTitleFont(42);
  eff.GetXaxis()->SetLabelFont(42);
  eff.GetYaxis()->SetTitleFont(42);
  eff.GetYaxis()->SetLabelFont(42);

  eff.GetXaxis()->SetLabelSize(0.055); // changed from 0.045
  eff.GetXaxis()->SetLabelOffset(0.008);
  eff.GetXaxis()->SetTickLength(0.03);
  eff.GetXaxis()->SetTitleSize(0.06); // changed from 0.05
  eff.GetXaxis()->SetTitleOffset(1.2);

  eff.GetYaxis()->SetTitleOffset(1.4); // changed from 1.8
  eff.GetYaxis()->SetTitleSize(0.06); // changed from 0.05
  eff.GetYaxis()->SetLabelSize(0.055); // changed from 0.045
  eff.GetYaxis()->SetTickLength(0.02);
  eff.GetYaxis()->SetLabelOffset(0.011);

  eff.GetYaxis()->SetRangeUser(0,1.5);
  eff.GetXaxis()->SetRangeUser(0,1000);

  eff.SetMarkerStyle(1);
  eff.SetLineWidth(2);
  eff.SetMarkerColor(1);


  eff.Draw("AP");

  TString infotext = TString::Format("%3.1f fb^{-1} (%d TeV)", 137., 13);
  TLatex *text = new TLatex(3.5, 24, infotext);
  text->SetNDC();
  text->SetTextAlign(33);
  text->SetX(0.95);
  text->SetTextFont(42);
  text->SetY(1);
  text->SetTextSize(0.045);
  text->Draw();

  m_can->SaveAs("plots/TrgEff_"+channel+".pdf");
}
