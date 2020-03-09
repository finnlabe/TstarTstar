// Efficiency histograms for trigger related studies
// author: F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q cutflowPlots.C

void cutflowPlot(){
    
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
  Double_t w = 600;
  Double_t h = 600;
  
  // Reuseable buffers
  TFile *input;
  TH1D *hist;

  // Drawing Definitions
  TCanvas *canvas = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.22,0.175,0.5,0.33);
  canvas->SetLogy();

  // Defining paths
  TString pathPresel = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Preselection/tgtg/hadded/";
  TString pathSel = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Selection/tgtg/hadded/";
  TString pathReco = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/MCStudy/tgtg/hadded/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";
  TString histname = "N_jets";

  // Defining Steps
  std::vector<TString> preselSteps = {"AfterCommon", "AfterMET"};
  std::vector<TString> selSteps = {"AfterMET", "After2D", "AfterNAK8sel_3"};
  std::vector<TString> recoSteps = {"AfterReco_Full"};
  int stepcount = preselSteps.size() + selSteps.size() + recoSteps.size();

  // Defining Samples
  std::vector<TString> signalSamples = {"TstarTstar_M-700.root", "TstarTstar_M-1000.root", "TstarTstar_M-1600.root"};
  std::vector<TString> BGSamples = {"TTbar.root", "WJets.root", "ST.root", "VV.root"};

  // Defining Drawing options
  std::vector<int> colors_Signal = {632, 820, 432};
  std::vector<int> colors_BG = {810, 800, 600, 416};
  std::vector<TString> labels = {"Initial", "Presel","MET", "2D", "N_AK8_Jet","AfterReco"};

  // ########################
  // ## Finish Definitions ##
  // ########################

  // Saving initial values
  std::vector<double> initial_signal;
  std::vector<double> initial_BG;
  double initial_BG_sum = 0;
  for(const auto & sample : signalSamples){
    input = TFile::Open(pathPresel+fileprefix+sample);
    hist = (TH1D*)input->Get(preselSteps.at(0)+"/"+histname);
    initial_signal.push_back(hist->Integral());
  }
  for(const auto & sample : BGSamples){
    input = TFile::Open(pathPresel+fileprefix+sample);
    hist = (TH1D*)input->Get(preselSteps.at(0)+"/"+histname);
    initial_BG.push_back(hist->Integral());
    initial_BG_sum += hist->Integral();
  }

  // Defining Hists
  std::vector<TH1*> cutflow_Signal;
  std::vector<TH1*> cutflow_BG;
  for(const auto & sample : signalSamples){cutflow_Signal.push_back(new TH1D(sample, sample, stepcount, 0, stepcount));}    // TODO remove ".root" in the end
  for(const auto & sample : BGSamples){cutflow_BG.push_back(new TH1D(sample, sample, stepcount, 0, stepcount));}
  
  // Filling hists
  int index_step = 1;

  // Preselection
  for(const auto & step : preselSteps){
    int index_sample = 0;
    for(const auto & sample : signalSamples){
      input = TFile::Open(pathPresel+fileprefix+sample);
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_Signal.at(index_sample)->SetBinContent(index_step, val/initial_signal.at(index_sample));
      index_sample++;
    }
    index_sample = 0;
    for(const auto & sample : BGSamples){
      input = TFile::Open(pathPresel+fileprefix+sample);
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_BG.at(index_sample)->SetBinContent(index_step, val/initial_BG_sum);
      index_sample++;
    }
    index_step++;
  }

  // Selection
  for(const auto & step : selSteps){
    int index_sample = 0;
    for(const auto & sample : signalSamples){
      input = TFile::Open(pathSel+fileprefix+sample);
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_Signal.at(index_sample)->SetBinContent(index_step, val/initial_signal.at(index_sample));
      index_sample++;
    }
    index_sample = 0;
    for(const auto & sample : BGSamples){
      input = TFile::Open(pathSel+fileprefix+sample);
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_BG.at(index_sample)->SetBinContent(index_step, val/initial_BG_sum);
      index_sample++;
    }
    index_step++;
  }

  // Reconstruction
  for(const auto & step : recoSteps){
    int index_sample = 0;
    for(const auto & sample : signalSamples){
      input = TFile::Open(pathReco+fileprefix+sample);
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_Signal.at(index_sample)->SetBinContent(index_step, val/initial_signal.at(index_sample));
      index_sample++;
    }
    index_sample = 0;
    for(const auto & sample : BGSamples){
      input = TFile::Open(pathReco+fileprefix+sample);
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_BG.at(index_sample)->SetBinContent(index_step, val/initial_BG_sum);
      index_sample++;
    }
    index_step++;
  }

  // ##########################
  // ## Finish Filling Hists ##
  // ##########################
  
  // drawing background
  THStack *BG_stack = new THStack("BG","");
  for(int index_draw_BG = cutflow_BG.size()-1; index_draw_BG >= 0; index_draw_BG--){
    int index_label = 1;
    for(const auto & label : labels) {cutflow_BG.at(index_draw_BG)->GetXaxis()->SetBinLabel(index_label, label); index_label++;}
    cutflow_BG.at(index_draw_BG)->SetFillColor(colors_BG.at(index_draw_BG)); 
    BG_stack->Add(cutflow_BG.at(index_draw_BG));
  }
  BG_stack->Draw("");

  // drawing signal
  int index_draw_signal = 0;
  for(const auto & hist : cutflow_Signal){
    hist->SetLineWidth(3); 
    hist->SetLineColor(colors_Signal.at(index_draw_signal)); // This can be improved by use of iterator
    hist->GetYaxis()->SetTitle("Efficiency");
    hist->Draw("same");
    index_draw_signal++;
  }
  
  // Draw line after certain bin to split presel, sel,
  bool doLine = false;
  if(doLine){
    double pos1 = cutflow_BG.at(0)->GetBinLowEdge(1+preselSteps.size());
    TLine *line1 = new TLine(pos1,1.25,pos1,1.9);
    line1->Draw("same");
    double pos2 = cutflow_BG.at(0)->GetBinLowEdge(1+preselSteps.size()+selSteps.size());
    TLine *line2 = new TLine(pos2,1.25,pos2,1.9);
    line2->Draw("same");
  }
  canvas->SaveAs("Cutflow.pdf");

  // Additional console oputput: efficiency of ttbarcut
  /**
  int index_sample = 0;
  for(const auto & sample : signalSamples){
    double pre_ttag = cutflow_Signal.at(index_sample)->GetBinContent(4);
    double post_ttag = cutflow_Signal.at(index_sample)->GetBinContent(5);
    cout << "For " << sample << " the efficiency of the t tag cut is " << post_ttag/pre_ttag << endl;
    index_sample++;
  }
  **/

}
