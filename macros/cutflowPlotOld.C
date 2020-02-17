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

  TString label="Cutflow";


  TString pathPresel = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Preselection/tgtg/hadded/";
  TString pathReco = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/MCStudy/tgtg/hadded/";
  TString fileprefix="uhh2.AnalysisModuleRunner.MC.";

  TString histname = "N_jets";

  int STEPCOUNT_Presel = 7;
  int STEPCOUNT_Reco = 1;
  TString stepsPresel[7] = {"AfterCommon", "AfterNjets", "AfterLepSel", "After2D", "AfterMET", "AfterST", "AfterAllCuts"};
  TString stepsReco[1] = {"AfterReco_Full"};
  
  int SIGNALSAMPLECOUNT = 3;
  int BGSAMPLECOUNT = 4;
  TString signalSamples[3] = {"TstarTstar_M-700.root", "TstarTstar_M-1000.root", "TstarTstar_M-1600.root"};
  TString BGSamples[4] = {"TTbar.root", "WJets.root", "ST.root", "VV.root"};

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.22,0.175,0.5,0.33);
  c1_hist->SetLogy();
  
  TFile *input;
  TH1D *hist;
  // Setting initial Integral values
  double signalInitial[3] = {0,0};
  for(int i = 0; i < SIGNALSAMPLECOUNT; i++){
    input = TFile::Open(pathPresel+fileprefix+signalSamples[i]);
    hist = (TH1D*)input->Get("AfterCommon/"+histname);
    signalInitial[i] = hist->Integral();
  }
  double BGInitial = 0;
  for(int i = 0; i < BGSAMPLECOUNT; i++){
    input = TFile::Open(pathPresel+fileprefix+BGSamples[i]);
    hist = (TH1D*)input->Get("AfterCommon/"+histname);
    BGInitial += hist->Integral();
  }
  
  TH1* cutflowBG = new TH1D("cutflowBG", "cutflowBG", STEPCOUNT_Presel+STEPCOUNT_Reco, 0, STEPCOUNT_Presel+STEPCOUNT_Reco);
  TH1* cutflowSignal[3] = {new TH1D("cutflow_M-700", "cutflow_M-700", STEPCOUNT_Presel+STEPCOUNT_Reco, 0, STEPCOUNT_Presel+STEPCOUNT_Reco),new TH1D("cutflow_M-1000", "cutflow_M-1000", STEPCOUNT_Presel+STEPCOUNT_Reco, 0, STEPCOUNT_Presel+STEPCOUNT_Reco), new TH1D("cutflow_M-1600", "cutflow_M-1600", STEPCOUNT_Presel+STEPCOUNT_Reco, 0, STEPCOUNT_Presel+STEPCOUNT_Reco)};

  // Presel
  for(int j = 0; j < STEPCOUNT_Presel; j++){
    // BG
    double BGval = 0;
    for(int i = 0; i < BGSAMPLECOUNT; i++){
      input = TFile::Open(pathPresel+fileprefix+BGSamples[i]);
      hist = (TH1D*)input->Get(stepsPresel[j]+"/"+histname);
      BGval += hist->Integral();
    }
    cutflowBG->GetXaxis()->SetBinLabel(j+1,stepsPresel[j]);
    cutflowBG->SetBinContent(j+1.5, BGval/BGInitial);

    // Signal
    double signalVal[3] = {0,0,0};
    for(int i = 0; i < SIGNALSAMPLECOUNT; i++){
      input = TFile::Open(pathPresel+fileprefix+signalSamples[i]);
      hist = (TH1D*)input->Get(stepsPresel[j]+"/"+histname);
      double Signalval = hist->Integral();
      cutflowSignal[i]->GetXaxis()->SetBinLabel(j+1,stepsPresel[j]);
      cutflowSignal[i]->SetBinContent(j+1, Signalval/signalInitial[i]);
    }
  }

  // Reco
  for(int j = 0; j < STEPCOUNT_Reco; j++){
    // BG
    double BGval = 0;
    for(int i = 0; i < BGSAMPLECOUNT; i++){
      input = TFile::Open(pathReco+fileprefix+BGSamples[i]);
      hist = (TH1D*)input->Get(stepsReco[j]+"/"+histname);
      BGval += hist->Integral();
    }
    cutflowBG->GetXaxis()->SetBinLabel(j+1+STEPCOUNT_Presel,stepsReco[j]);
    cutflowBG->SetBinContent(j+1+STEPCOUNT_Presel, BGval/BGInitial);

    // Signal
    double signalVal[3] = {0,0,0};
    for(int i = 0; i < SIGNALSAMPLECOUNT; i++){
      input = TFile::Open(pathReco+fileprefix+signalSamples[i]);
      hist = (TH1D*)input->Get(stepsReco[j]+"/"+histname);
      double Signalval = hist->Integral();
      cutflowSignal[i]->GetXaxis()->SetBinLabel(j+1+STEPCOUNT_Presel,stepsReco[j]);
      cutflowSignal[i]->SetBinContent(j+1+STEPCOUNT_Presel, Signalval/signalInitial[i]);
    }
  }

  leg->AddEntry(cutflowBG,"Background","l");
  cutflowBG->SetTitle("Cutflow");
  cutflowBG->SetLineWidth(3); 
  cutflowBG->SetLineColor(1);
  cutflowBG->GetYaxis()->SetTitle("Efficiency");
  cutflowBG->Draw("");
  int signalColors[3] = {632, 416, 432};
  for(int i = 0; i < SIGNALSAMPLECOUNT; i++){
    leg->AddEntry(cutflowSignal[i],signalSamples[i],"l");
    cutflowSignal[i]->SetLineWidth(3); 
    cutflowSignal[i]->SetLineColor(signalColors[i]);   
    cutflowSignal[i]->GetYaxis()->SetTitle("Efficiency");
    cutflowSignal[i]->Draw("same");
  }

  leg->Draw("same");
  c1_hist->SaveAs(label+".pdf");

}
