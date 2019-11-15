// Efficiency histograms for trigger related studies
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q DeltaRPlots.C

void DeltaRPlots(TString filename="uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgammaTgluon_M-700_Run2016v3.root", TString label="TstarTstarToTgammaTgluon_", TString subpath="RecoPlots_After_ttbar",TString histname="DeltaR_tophad_gamma"){
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.06,"x");  
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"x");  
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetPalette(55);

  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.17);

  gROOT->ForceStyle();
  Double_t w = 1200;
  Double_t h = 1000;

  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Preselection/RunII_2016_MuonHihjPtId_mvaPhoIDwp90Fall17_nonIsoandIsoHLT_addTTBarRECO/semileptonic/";
  TString filenamepre = "uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgammaTgluon_M-";
  TString filenamepost = "_Run2016v3.root";

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.5,0.7,0.83,0.90);

  TString masspoints[] = {"700", "1000", "1300", "1600"};

  int i = 2;
  for(auto masspoint: masspoints){
    TFile *input = TFile::Open(path+filenamepre+masspoint+filenamepost);
    TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname);//histogram
    if(!hist) cout<<"Hist is empty"<<endl;;
    hist->GetYaxis()->SetRangeUser(0., 0.02);
    hist->GetXaxis()->SetTitle("#DeltaR");
    hist->GetYaxis()->SetTitle("events");
    hist->SetTitle(histname);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(i);
    hist->SetLineColor(i);
    hist->Draw("hist same");

    leg->AddEntry(hist,"M_{T*} = "+masspoint,"l");

    i++;
  }
  
  leg->Draw("same");
  c1_hist->SaveAs("DeltaR_"+label+"_"+subpath+"_"+histname+".pdf");
}
