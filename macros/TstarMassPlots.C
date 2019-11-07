// Efficiency histograms for trigger related studies
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q photonSelecPlots.C

void TstarMassPlots(TString filename="uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgammaTgluon_M-1500_Run2016v3.root", TString label="TstarTstarToTgammaTgluon_M-1500_Run2016v3", TString subpath="AfterMatching",TString histname="M_Tstar_gluon"){
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
  filename = "uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgammaTgluon_M-1500_Run2016v3";
  
  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.5,0.7,0.83,0.95);

  TFile *input = TFile::Open(path+"/old/"+filename+"___noCut.root");
  TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname);//histogram
  if(!hist) cout<<"Hist is empty"<<endl;;
  hist->GetXaxis()->SetTitle("M_{T*} (reconstructed)");
  hist->GetYaxis()->SetTitle("events");
  hist->SetTitle("M_{T^{*}} #DeltaR Cut Comparison");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(1);
  hist->SetLineColor(2);
  hist->Draw("hist");
  leg->AddEntry(hist,"no cut","pl");    


  TFile *input2 = TFile::Open(path+"/old/"+filename+"___02cut.root");
  TH1D *hist2 = (TH1D*)input->Get(subpath+"/"+histname);//histogram
  if(!hist2) cout<<"Hist is empty"<<endl;;
  hist2->GetXaxis()->SetTitle("M^{gluonic}_{T*} (reconstructed)");
  hist2->GetYaxis()->SetTitle("events");
  hist2->SetTitle("M_{T^{*}} #DeltaR Cut Comparison");
  hist2->SetMarkerStyle(20);
  hist2->SetMarkerColor(2);
  hist2->SetLineColor(2);
  hist2->Draw("hist");
  leg->AddEntry(hist2,"0.2 cut","pl");    
  
  TFile *input3 = TFile::Open(path+"/old/"+filename+"___1cut.root");
  TH1D *hist3 = (TH1D*)input->Get(subpath+"/"+histname);//histogram
  if(!hist3) cout<<"Hist is empty"<<endl;;
  hist3->GetXaxis()->SetTitle("M^{gluonic}_{T*} (reconstructed)");
  hist3->GetYaxis()->SetTitle("events");
  hist3->SetTitle("M_{T^{*}} #DeltaR Cut Comparison");
  hist3->SetMarkerStyle(20);
  hist3->SetMarkerColor(3);
  hist3->SetLineColor(3);
  hist3->Draw("hist");
  leg->AddEntry(hist3,"1 cut","pl");    

  leg->Draw();
 
  c1_hist->SaveAs("MTstar_"+label+"_"+subpath+"_"+histname+"_DeltaRCutcomparison.pdf");
}
