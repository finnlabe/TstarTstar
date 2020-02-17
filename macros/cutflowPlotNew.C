// Efficiency histograms for trigger related studies
// author: F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q cutflowPlots.C

void cutflowPlotNew(){
    
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

  TString label="Cutflow";


  TString pathPresel = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Preselection/tgtg/hadded/";
  TString pathSel = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Selection/tgtg/hadded/";
  TString pathReco = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/MCStudy/tgtg/hadded/";
  TString fileprefix="uhh2.AnalysisModuleRunner.MC.";

  TString histname = "N_jets";

  TString stepsPresel[7] = {"AfterCommon", "AfterNjets", "AfterLepSel", "After2D", "AfterMET", "AfterST", "AfterAllCuts"};
  TString stepsSel[0];
  TString stepsReco[1] = {"AfterReco_Full"};
  
  TString signalSamples[3] = {"TstarTstar_M-700.root", "TstarTstar_M-1000.root", "TstarTstar_M-1600.root"};
  TString BGSamples[4] = {"TTbar.root", "WJets.root", "ST.root", "VV.root"};

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.22,0.175,0.5,0.33);
  c1_hist->SetLogy();
  
  TFile *input;
  TH1D *hist;

  // Setting initial Integral values
  double signalInitial[3] = {0,0,0};
  for(int i = 0; i < 3; i++){
    input = TFile::Open(pathPresel+fileprefix+signalSamples[i]);
    hist = (TH1D*)input->Get("AfterCommon/"+histname);
    signalInitial[i] = hist->Integral();
  }
  double BGInitial = 0;
  for(int i = 0; i < 4; i++){
    input = TFile::Open(pathPresel+fileprefix+BGSamples[i]);
    hist = (TH1D*)input->Get("AfterCommon/"+histname);
    BGInitial += hist->Integral();
  }
  
  TH1* cutflowBG[4] = {new TH1D("cutflowTTBar", "cutflowTTBar", 8, 0, 8), new TH1D("cutflowWJets", "cutflowWJets", 8, 0, 8), new TH1D("cutflowST", "cutflowST", 8, 0, 8),new TH1D("cutflowVV", "cutflowVV", 8, 0, 8)};
  TH1* cutflowSignal[3] = {new TH1D("cutflow_M-700", "cutflow_M-700", 8, 0, 8),new TH1D("cutflow_M-1000", "cutflow_M-1000", 8, 0, 8), new TH1D("cutflow_M-1600", "cutflow_M-1600", 8, 0, 8)};

  // Presel
  for(int j = 0; j < 7; j++){
    // BG
    for(int i = 0; i < 4; i++){
      input = TFile::Open(pathPresel+fileprefix+BGSamples[i]);
      hist = (TH1D*)input->Get(stepsPresel[j]+"/"+histname);
      double BGval = hist->Integral();
      cutflowBG[i]->GetXaxis()->SetBinLabel(j+1,stepsPresel[j]);
      cutflowBG[i]->SetBinContent(j+1.5, BGval/BGInitial);
    }

    // Signal
    for(int i = 0; i < 3; i++){
      input = TFile::Open(pathPresel+fileprefix+signalSamples[i]);
      hist = (TH1D*)input->Get(stepsPresel[j]+"/"+histname);
      double Signalval = hist->Integral();
      cutflowSignal[i]->GetXaxis()->SetBinLabel(j+1,stepsPresel[j]);
      cutflowSignal[i]->SetBinContent(j+1, Signalval/signalInitial[i]);
    }
  }

  // Sel

  // Reco
  for(int j = 0; j < 1; j++){
    // BG
    for(int i = 0; i < 4; i++){
      input = TFile::Open(pathReco+fileprefix+BGSamples[i]);
      hist = (TH1D*)input->Get(stepsReco[j]+"/"+histname);
      double BGval = hist->Integral();
      cutflowBG[i]->GetXaxis()->SetBinLabel(j+1+7,stepsReco[j]);
      cutflowBG[i]->SetBinContent(j+1+7, BGval/BGInitial);
    }
    
    // Signal
    for(int i = 0; i < 3; i++){
      input = TFile::Open(pathReco+fileprefix+signalSamples[i]);
      hist = (TH1D*)input->Get(stepsReco[j]+"/"+histname);
      double Signalval = hist->Integral();
      cutflowSignal[i]->GetXaxis()->SetBinLabel(j+1+7,stepsReco[j]);
      cutflowSignal[i]->SetBinContent(j+1+7, Signalval/signalInitial[i]);
    }
  }
  
  int BGColors[4] = {810, 800, 600, 416};
  int signalColors[3] = {632, 416, 432};

  THStack *stackBG = new THStack("BG","");
  for(int i = 3; i >= 0; i--){
    cutflowBG[i]->SetFillColor(BGColors[i]);
    stackBG->Add(cutflowBG[i]);
  }

  stackBG->Draw();

  for(int i = 0; i < 3; i++){
    leg->AddEntry(cutflowSignal[i],signalSamples[i],"l");
    cutflowSignal[i]->SetLineWidth(3); 
    cutflowSignal[i]->SetLineColor(signalColors[i]);   
    cutflowSignal[i]->GetYaxis()->SetTitle("Efficiency");
    cutflowSignal[i]->Draw("same");
  }

  leg->Draw("same");
  c1_hist->SaveAs(label+".pdf");

}
