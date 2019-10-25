// Efficiency histograms for trigger related studies
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q photonSelecPlots.C

void TstarMassPlots(TString filename="uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgammaTgluon_M-1500_Run2016v3.root", TString label="TstarTstarToTgammaTgluon_M-1500_Run2016v3", TString subpath="AfterMatching",TString histname="M_Tstar"){
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.06,"x");  
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"x");  
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetPalette(55);

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.17);

  gROOT->ForceStyle();
  Double_t w = 1200;
  Double_t h = 600;

  //TODO: should use some kind of for loop if this is actually used in any way!
  //Files after selection
  
  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Preselection/RunII_2016_MuonHihjPtId_mvaPhoIDwp90Fall17_nonIsoandIsoHLT_addTTBarRECO/semileptonic/";
  
  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  
  TFile *input = TFile::Open(path+filename);
  TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname);//histogram
  if(!hist) cout<<"Hist is empty"<<endl;;
  hist->GetXaxis()->SetTitle("M_{T*} (reconstructed)");
  hist->GetYaxis()->SetTitle("events");
  hist->GetYaxis()->SetRangeUser(0,400);
  hist->SetTitle("");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(1);
  hist->Draw();
 
  c1_hist->SaveAs("MTstar_"+label+"_"+subpath+"_"+histname+".pdf");
}
