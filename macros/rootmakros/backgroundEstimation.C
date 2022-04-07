// author: F.Labe
// date: 24.09.2021
// Run it with following command:
// root -l -b -q backgroundEstimation.C

template<class T, size_t N>
constexpr size_t size(T (&)[N]) { return N; }

// stolen from alex :)
void fill_reweighted(TH1F* hReweight,TH1F* hReweight_up,TH1F* hReweight_down, TF1* fitfun1, TF1* fitfun2,  TFitResultPtr r1, TFitResultPtr r2, const TString &fname, const double &weightSign) {

  TFile *f = TFile::Open(fname);
  if(!f) std::cout << "File not found!" << std::endl;

  TTreeReader tReader("AnalysisTree", f);
  TTreeReaderValue<double> preWeight(tReader, "evt_weight");
  TTreeReaderValue<double> recoST(tReader, "ST");
  TTreeReaderValue<TString> regionReader(tReader, "region");

  while (tReader.Next()) {
    double ST = (*recoST);
    TString region = (*regionReader);
    if(region != "CR2") continue;

    double x[1] = { ST };
    double err_1[1], err_2[1];
    // mass = (mass >= 4000.) ? 3999.9 : mass; // I think I dont need an overflow bin after rebinning, but maybe check!

    // evaluate fit at masspoint and get 1 sigma deviation from 68.3% confidence interval
    double diff = fitfun1->Eval(*recoST) - fitfun2->Eval(*recoST);
    double weight = fitfun1->Eval(*recoST);
    r1->GetConfidenceIntervals(1, 1, 1, x, err_1, 0.683, false);
    r2->GetConfidenceIntervals(1, 1, 1, x, err_2, 0.683, false);
    double err = sqrt( pow(diff/2, 2) + pow(err_1[0], 2) + pow(err_2[0], 2));
    double weight_up = weight - diff/2 + err;
    double weight_down = weight - diff/2 - err;
    weight = (weight < 0) ? 0 : weight * (*preWeight);
    weight_up = (weight_up < 0) ? 0 : weight_up * (*preWeight);
    weight_down = (weight_down < 0) ? 0 : weight_down * (*preWeight);
    hReweight->Fill(ST, weight * weightSign);
    hReweight_up->Fill(ST, weight_up * weightSign);
    hReweight_down->Fill(ST, weight_down * weightSign);
  }

  f->Close();
  delete f;
}

// stolen also from alex :)
void create_output(const TString fout_name, const TString subdir_name, TH1F* hist) {
  TFile *fout = new TFile(fout_name, "update");
  TDirectory *subdir = (TDirectory*) fout->Get(subdir_name); if (!subdir) subdir = fout->mkdir(subdir_name);
  // TDirectory* subdir = fout->mkdir(subdir_name);
  subdir->cd();
  TString hname = "pt_ST_rebinned";
  TH1F *hOut =(TH1F*) hist->Clone(hname);
  hOut->Write();
  fout->Close();
  delete fout;
}


void backgroundEstimation(){

  // definitions
  std::vector<TString> nontop_backgrounds = {"WJets", "QCD", "VV"};
  std::vector<TString> top_backgrounds = {"ST", "TTbar"};
  TString data = "data";

  TString subpath_SR="newTaggerCR";
  TString subpath_CR="newTagger_btagCR";
  TString histname="pt_ST_rebinned";
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/UL18/hadded/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.";

  TH1D *histSR;
  TH1D *hist_btagCR_nontop;
  TH1D *hist_btagCR_top;
  TH1D *hist_btagCR;
  TH1D *hist_btagCR_data;

  // open full set of non-top backgrounds in signal region
  bool first = true;
  for (auto filename : nontop_backgrounds) {
    TFile *input = TFile::Open(path+fileprefix+"MC."+filename+".root");
    if(!input) cout << "Empty file" << endl;
    TH1D *hist = (TH1D*)input->Get(subpath_SR+"/"+histname);
    if(!hist) cout << "Empty hist" << endl;
    if (first) {
      first = false;
      histSR = hist;
    } else {
      histSR->Add(hist);
    }
  }

  // open full set of non-top backgrounds in W-jets control region
  first = true;
  vector<double> nontop_integrals;
  for (auto filename : nontop_backgrounds) {
    TFile *input = TFile::Open(path+fileprefix+"MC."+filename+".root");
    if(!input) cout << "Empty file" << endl;
    TH1D *hist = (TH1D*)input->Get(subpath_CR+"/"+histname);
    if(!hist) cout << "Empty hist" << endl;
    // outputting some information on the b-tagging CR
    nontop_integrals.push_back(hist->Integral());
    if (first) {
      first = false;
      hist_btagCR_nontop = hist;
    } else {
      hist_btagCR_nontop->Add(hist);
    }
  }

  // open the top backgrounds
  first = true;
  vector<double> top_integrals;
  for (auto filename : top_backgrounds) {
    TFile *input = TFile::Open(path+fileprefix+"MC."+filename+".root");
    if(!input) cout << "Empty file" << endl;
    TH1D *hist = (TH1D*)input->Get(subpath_CR+"/"+histname);
    if(!hist) cout << "Empty hist" << endl;
    // outputting some information on the b-tagging CR
    top_integrals.push_back(hist->Integral());
    if (first) {
      first = false;
      hist_btagCR_top = hist;
    } else {
      hist_btagCR_top->Add(hist);
    }
  }

  // combining for full picture
  hist_btagCR = (TH1D*) hist_btagCR_nontop->Clone();
  hist_btagCR->Add(hist_btagCR_top);

  for(uint i = 0; i < nontop_backgrounds.size(); i++) std::cout << "Contribution of " << nontop_backgrounds.at(i) << ": " << nontop_integrals.at(i)/hist_btagCR->Integral() << std::endl;
  for(uint i = 0; i < top_backgrounds.size(); i++) std::cout << "Contribution of " << top_backgrounds.at(i) << ": " << top_integrals.at(i)/hist_btagCR->Integral() << std::endl;

  TGraphAsymmErrors purity = TGraphAsymmErrors();
  purity.Divide(hist_btagCR_nontop, hist_btagCR, "cl=0.68 b(1,1) mode");

  TFile *file = new TFile("files/bgest_purity.root", "RECREATE");
  purity.SetName("purity");
  purity.Write();
  file->Close();
  delete file;

  // calculate ratio histogram
  TGraphAsymmErrors ratio = TGraphAsymmErrors();
  ratio.Divide(histSR, hist_btagCR, "pois");

  // fit function to ratio histogram
  TF1 *fit = new TF1("fit", "pol 1", 500, 6000);
  TH1D *hint1 = new TH1D("hint", "Fit 1 with conf.band", 100, 0, 10000);
  TF1 *fit2 = new TF1("fit2", "pol 2", 500, 6000);
  TH1D *hint2 = new TH1D("hint", "Fit 2 with conf.band", 100, 0, 10000);
  TFitResultPtr r_1 = ratio.Fit("fit", "NS", "", 500, 6000);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint1, 0.68);
  TFitResultPtr r_2 = ratio.Fit("fit2", "NS", "", 500, 6000);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint2, 0.68);

  const int nbins = 34;
  double bins[nbins] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500,
    2600, 2700, 2800, 2900, 3000, 3250, 4000, 6000};
  TH1F *hReweight_fit_M = new TH1F("reweight_fit", "_reweight_fit", nbins-1, bins);
  TH1F *hReweight_fit_M_up = new TH1F("reweight_fit_up", "_reweight_fit", nbins-1, bins);
  TH1F *hReweight_fit_M_down = new TH1F("reweight_fit_down", "_reweight_fit", nbins-1, bins);
  hReweight_fit_M->Sumw2();
  hReweight_fit_M_up->Sumw2();
  hReweight_fit_M_down->Sumw2();
  fill_reweighted(hReweight_fit_M, hReweight_fit_M_up, hReweight_fit_M_down, fit, fit2, r_1, r_2, path + fileprefix+ "DATA.DATA.root", 1);

  for (int i = 0; i < top_backgrounds.size(); ++i){
    TString fName = path + fileprefix + "MC." + top_backgrounds.at(i) + ".root";
    fill_reweighted(hReweight_fit_M, hReweight_fit_M_up, hReweight_fit_M_down, fit, fit2, r_1, r_2, fName, -1);
  }

  const TString prefix = "uhh2.AnalysisModuleRunner.";
  create_output(prefix+"MC.Other.root", "reco", hReweight_fit_M);
  create_output(prefix+"MC.Other_UP.root", "reco", hReweight_fit_M_up);
  create_output(prefix+"MC.Other_DOWN.root", "reco", hReweight_fit_M_down);

  // plot everything
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

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  c1_hist->SetLogz();

  // draw plot
  ratio.GetXaxis()->SetTitle("S_{T} [GeV]");
  ratio.GetXaxis()->SetNdivisions(505);
  ratio.GetYaxis()->SetTitle("ratio");
  ratio.GetYaxis()->SetRangeUser(1, 35);
  ratio.SetTitle("");

  ratio.Draw("");
  fit->Draw("same");
  fit2->SetLineColor(4);
  fit2->Draw("same");

  hint1->SetLineColor(2);
  hint1->SetFillColorAlpha(2,0.3);
  hint1->Draw("e3 same");
  hint2->SetFillColorAlpha(4,0.3);
  hint1->SetLineColor(4);
  hint2->Draw("e3 same");

  // draw Lumi text
  TString infotext = TString::Format("%3.0f fb^{-1} (%d TeV)", 138., 13);
  TLatex *text = new TLatex(3.5, 24, infotext);
  text->SetNDC();
  text->SetTextAlign(33);
  text->SetX(0.88);
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
  text3->SetX(0.29);
  text3->SetTextFont(52);
  text3->SetTextSize(0.035);
  text3->SetY(0.986);
  text3->Draw();

  c1_hist->SaveAs("plots/backgroundEstimation.pdf");
  c1_hist->Clear();

  purity.GetXaxis()->SetTitle("S_{T} [GeV]");
  purity.GetXaxis()->SetNdivisions(505);
  purity.GetYaxis()->SetTitle("purity");
  purity.SetTitle("");

  purity.Draw("AP");

  text->Draw();
  text2->Draw();
  text3->Draw();

  c1_hist->SaveAs("plots/purity.pdf");


}
