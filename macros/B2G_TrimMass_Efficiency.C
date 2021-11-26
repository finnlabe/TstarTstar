// this makro takes the "before" hist (from before_path) and the "after Trigger" hists (afterPaths) and calculates efficiencies by deviding.
// Everything else is just input reading and a lot of styling options

void B2G_TrimMass_Efficiency(){

  gROOT->SetBatch(kTRUE); // to not open canvas and get XQuartz in the way

  // defining histogram to read
  TString histname = "pt";
  bool doSignalSample = false;
  bool doComparison = false; // only one can be true, doComparison is stronger

  // input paths, filenames...
  TString path = "/nfs/dust/cms/user/flabe/B2G_Trigger_contact/data/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.";
  TString filename = "DATA.SingleMuon2018_RunD_TrimMassTriggerTest_2018";
  TString signalSampleFilename = "MC.NMSSM_UL18";
  TString before_path = "TrimMass_before";
  std::vector<TString> afterPaths = {"TrimMass_afterHT", "TrimMass_afterPFJet500", "TrimMass_TrimMass30", "TrimMass_TrimMass50", "TrimMass_afterAll"};
  std::vector<TString> titles = {"PFHT_1050", "PFJet500", "AK8PFJet400_TrimMass30", "AK8PFHT800_TrimMass50", "Total"};
  std::vector<int> colors = {kGreen+2,kRed,kBlue,kBlue-9,1};

  if(doComparison) {
    doSignalSample = false;
    afterPaths = {"TrimMass_afterAll", "TrimMass_currentPlan", "TrimMass_PFJet420_TrimMass30", "TrimMass_PFJet500_TrimMass30"};
    titles = {"2018 setup", "Current plan", "with PFJet420_TrimMass30", "with PFJet500_TrimMass30"};
    colors = {1, kOrange+7, kViolet, kViolet-6};
  }

  // set x label
  TString xAxisLabel;
  if(histname.BeginsWith("HT")) xAxisLabel = "H_{T} [GeV]";
  if(histname.BeginsWith("mjj")) xAxisLabel = "m_{jj} [GeV]";
  else if(histname.BeginsWith("mj")) xAxisLabel = "m_{j} [GeV]";
  if(histname.BeginsWith("pt")) xAxisLabel = "p_{T} [GeV]";

  // open file
  TFile *file = TFile::Open(path+fileprefix+filename+".root");

  // doing the quick alpha check
  /**
  TH1D *hist_alpha_ref = (TH1D*)file->Get("h_alpha_ref/oneBin");
  TH1D *hist_alpha_mu = (TH1D*)file->Get("h_alpha_mu/oneBin");
  TH1D *hist_alpha_jets = (TH1D*)file->Get("h_alpha_jets/oneBin");
  TH1D *hist_alpha_both= (TH1D*)file->Get("h_alpha_both/oneBin");

  double eff_mu = hist_alpha_mu->Integral() / hist_alpha_ref->Integral();
  double eff_jets = hist_alpha_jets->Integral() / hist_alpha_ref->Integral();
  double eff_both = hist_alpha_both->Integral() / hist_alpha_ref->Integral();
  std::cout << "Efficiencies: " << std::endl;
  std::cout << "Single muon trigger (reference): " << eff_mu << std::endl;
  std::cout << "Combined jets triggers: " << eff_jets << std::endl;
  std::cout << "Both together: " << eff_both << std::endl;

  double alpha = (eff_mu * eff_jets)/eff_both;
  std::cout << "The alpha parameter is " << alpha << "." << std::endl;
  **/

  // read before hist histograms
  TH1D *hist_before = (TH1D*)file->Get(before_path+"/"+histname);

  // plotting parameters
  Double_t w = 600;
  Double_t h = 400;
  TCanvas *canvas = new TCanvas("chist", "c", w, h);

  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.16);  gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.09);
  gPad->SetTicky();

  TLegend *leg = new TLegend(0.3,0.2);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);

  // calculate and plot efficiency hists
  std::vector<TGraphAsymmErrors> effs = {TGraphAsymmErrors(), TGraphAsymmErrors(), TGraphAsymmErrors(), TGraphAsymmErrors(), TGraphAsymmErrors()};
  int i = 0;
  for (const auto & subpath : afterPaths) {
    std::cout << "Processing " << subpath << "." << std::endl;
    TH1D *hist = (TH1D*)file->Get(subpath+"/"+histname);

    // this line is the actual efficiency calculation!
    effs.at(i).Divide(hist, hist_before, "cl=0.68 b(1,1) mode");

    effs.at(i).GetXaxis()->SetTitle(xAxisLabel);
    effs.at(i).GetYaxis()->SetTitle("Efficiency");
    effs.at(i).SetTitle("");
    effs.at(i).SetMarkerStyle(1);
    effs.at(i).SetLineWidth(2);
    effs.at(i).SetLineColor(colors.at(i));
    effs.at(i).SetMarkerColor(colors.at(i));

    // Drawing
    if(i == 0) {
      // cosmetics
      effs.at(i).GetXaxis()->SetTitleFont(42);
      effs.at(i).GetXaxis()->SetLabelFont(42);
      effs.at(i).GetYaxis()->SetTitleFont(42);
      effs.at(i).GetYaxis()->SetLabelFont(42);

      effs.at(i).GetXaxis()->SetLabelSize(0.055); // changed from 0.045
      effs.at(i).GetXaxis()->SetLabelOffset(0.008);
      effs.at(i).GetXaxis()->SetTickLength(0.03);
      effs.at(i).GetXaxis()->SetTitleSize(0.06); // changed from 0.05
      effs.at(i).GetXaxis()->SetTitleOffset(1.2);

      effs.at(i).GetYaxis()->SetTitleOffset(0.9); // changed from 1.8
      effs.at(i).GetYaxis()->SetTitleSize(0.06); // changed from 0.05
      effs.at(i).GetYaxis()->SetLabelSize(0.055); // changed from 0.045
      effs.at(i).GetYaxis()->SetTickLength(0.02);
      effs.at(i).GetYaxis()->SetLabelOffset(0.011);

      effs.at(i).GetXaxis()->SetNdivisions(505);
      effs.at(i).GetYaxis()->SetRangeUser(0,1.5);
      //effs.at(i).GetXaxis()->SetRangeUser(0,2000);

      effs.at(i).Draw("ap");
    }
    else {
      effs.at(i).Draw("p same");
    }
    leg->AddEntry(&effs.at(i), titles.at(i), "l");

    std::cout << "Done." << std::endl;
    i++;
  }

  leg->SetBorderSize(0);
  leg->Draw("same");

  // draw Lumi text
  TString infotext = TString::Format("Run 2018D (%d TeV)", 13);
  TLatex *text = new TLatex(3.5, 24, infotext);
  text->SetNDC();
  text->SetTextAlign(33);
  text->SetX(0.91);
  text->SetTextFont(42);
  text->SetY(1);
  text->SetTextSize(0.045);
  text->Draw();

  // draw CMS Work in Progress text
  TString cmstext = "CMS";
  TLatex *text2 = new TLatex(3.5, 24, cmstext);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.185);
  text2->SetTextFont(62);
  text2->SetTextSize(0.05);
  text2->SetY(1);
  text2->Draw();
  TString preltext = "Work in Progress";
  TLatex *text3 = new TLatex(3.5, 24, preltext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.25);
  text3->SetTextFont(52);
  text3->SetTextSize(0.035);
  text3->SetY(0.99);
  text3->Draw();

  // write cuts in plots additional text
  // TODO finish this
  double position = 0.4;
  double step = 0.05;
  double posX = 0.7;
  double size = 0.035;

  TString cuttext1;
  if(histname.Contains("600pt")) cuttext1 = "p_{T}^{AK8} > 600 GeV";
  else cuttext1 = "p_{T}^{AK8} > 200 GeV";
  TLatex *textCuts1 = new TLatex(3.5, 24, cuttext1);
  textCuts1->SetNDC();
  textCuts1->SetX(posX);
  textCuts1->SetTextSize(size);
  textCuts1->SetY(position);
  textCuts1->Draw();
  position-=step;

  TString cuttext2 = "|#eta^{AK8}| < 2.4";
  TLatex *textCuts2 = new TLatex(3.5, 24, cuttext2);
  textCuts2->SetNDC();
  textCuts2->SetX(posX);
  textCuts2->SetTextSize(size);
  textCuts2->SetY(position);
  textCuts2->Draw();
  position-=step;

  TString cuttext3 = "#Delta#eta_{jj} < 1.3";
  TLatex *textCuts3 = new TLatex(3.5, 24, cuttext3);
  textCuts3->SetNDC();
  textCuts3->SetX(posX);
  textCuts3->SetTextSize(size);
  textCuts3->SetY(position);
  textCuts3->Draw();
  position-=step;

  if(histname.Contains("abovemj")) {
    TString cuttext4;
    if(histname.Contains("55")) cuttext4 = "m_{j} > 55 GeV";
    else cuttext4 = "m_{j} > 35 GeV";
    TLatex *textCuts4 = new TLatex(3.5, 24, cuttext4);
    textCuts4->SetNDC();
    textCuts4->SetX(posX);
    textCuts4->SetTextSize(size);
    textCuts4->SetY(position);
    textCuts4->Draw();
    position-=step;
  }

  if(histname.Contains("aboveHT")) {
    TString cuttext5 = "H_{T} > 1000 GeV";
    TLatex *textCuts5 = new TLatex(3.5, 24, cuttext5);
    textCuts5->SetNDC();
    textCuts5->SetX(posX);
    textCuts5->SetTextSize(size);
    textCuts5->SetY(position);
    textCuts5->Draw();
    position-=step;
  }

  // plotting a signal sample in here
  if (doSignalSample) {
    TFile *signalfile = TFile::Open(path+fileprefix+signalSampleFilename+".root");
    TObjArray *tx = histname.Tokenize("_");
    TString signalSampleHistname = ((TObjString *)(tx->At(0)))->String();
    std::cout << signalSampleHistname << std::endl;
    TH1D *signalSampleHist = (TH1D*)signalfile->Get("main/"+signalSampleHistname);
    if(!signalSampleHist) std::cout << "Hist does not exist" << std::endl;

    // scale signalSampleHist to the pad coordinates
    signalSampleHist->Scale(1/signalSampleHist->Integral());
    signalSampleHist->SetFillColorAlpha(14, 0.5);
    signalSampleHist->SetLineColor(0);
    signalSampleHist->Draw("hist same");

  }

  TString outputname = "plots/B2G_TrimMass_Eff_"+histname+".pdf";
  if(doSignalSample) outputname = "plots/B2G_TrimMass_eff_"+histname+"_withSignal.pdf";
  if(doComparison) outputname = "plots/B2G_TrimMass_comparison_"+histname+".pdf";
  canvas->SaveAs(outputname);

}
