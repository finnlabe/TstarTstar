// Just a basic macro to create ratio plot of two hists.
// author: F. Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q ratioplotter.C

void ULcompare_ratioplotter(){

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
  Double_t w = 800;
  Double_t h = 600;

  TString sample = "TTbar";

  TString EOYpath = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/2018/hadded/";
  TString ULpath = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/UL18/hadded/";
  TString filename = "uhh2.AnalysisModuleRunner.MC."+sample+".root";
  TString histname = "pt_ST";
  TString label = "S_{T}";

  TString folder = "AfterST";

  TCanvas *canvas = new TCanvas("canvas", "c", w, h);
  canvas->SetLogy();

  TFile *input_EOY = TFile::Open(EOYpath+filename);
  TFile *input_UL = TFile::Open(ULpath+filename);

  TH1D *hist_num = (TH1D*)input_UL->Get(folder+"/"+histname); //histogram
  TH1D *hist_denom = (TH1D*)input_EOY->Get(folder+"/"+histname); //histogram
  if(!hist_num) std::cout << "Numerator does not exist!" << std::endl;
  if(!hist_denom) std::cout << "Denominator does not exist!" << std::endl;

  hist_num->SetTitle("");
  hist_num->GetXaxis()->SetTitle(label);
  hist_num->GetYaxis()->SetTitle("events");
  hist_num->SetMarkerStyle(20);
  hist_num->SetMarkerColor(1);
  hist_num->SetLineColor(1);
  hist_num->Draw("hist");
  hist_denom->SetMarkerStyle(20);
  hist_denom->SetMarkerColor(2);
  hist_denom->SetLineColor(2);
  hist_denom->Draw("hist same");

  auto legend = new TLegend(0.55,0.7,0.78,0.88);
  gStyle->SetLegendTextSize(0.05);
  legend->AddEntry(hist_num,"UL","l");
  legend->AddEntry(hist_denom,"EOY","l");
  legend->Draw();

  canvas->SaveAs("plots/ULcomparison_"+histname+"_"+sample+".pdf");

  canvas->SetLogy(false);

  // plot ratio
  hist_num->Divide(hist_denom);
  hist_num->GetYaxis()->SetTitle("ratio");
  hist_num->SetLineColor(2);
  hist_num->Draw("hist");
  canvas->SaveAs("plots/STratio_"+histname+"_"+sample+".pdf");

}
