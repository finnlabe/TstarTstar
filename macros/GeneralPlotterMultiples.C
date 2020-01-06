// Just a basic macro to create simple plots from given hists.
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q GeneralPlotterMultiples.C

void GeneralPlotterMultiples(TString channel = "tgtgamma", TString masspoint="1600", TString histname="Mtophad_bestchi2"){
  
  const int count = 2;
  TString subpaths[count] = {"RecoPlots_GEN", "RecoPlots_After_TstarTstar"};

  Double_t w = 800;
  Double_t h = 600;

  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Preselection/RunII_2016_MuonHihjPtId_mvaPhoIDwp90Fall17_nonIsoandIsoHLT_addTTBarRECO/"+channel+"/old/GEN_first/";
  TString filename = "uhh2.AnalysisModuleRunner.MC.MC_TstarTstarTo";
  if(channel == "tgtg"){filename += "TgluonTgluon";}
  else if(channel == "tgtgamma"){filename += "TgammaTgluon";}
  else {cout<<"Channel invalid!"<<endl;}
  filename += "_M-"+masspoint+"_Run2016v3.root";
  
  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);

  for(int i = 0; i < count; i++){
    TFile *input = TFile::Open(path+filename);
    TH1D *hist = (TH1D*)input->Get(subpaths[i]+"/"+histname); //histogram
    hist->Scale(1/hist->Integral());
    if(!hist) cout<<"Hist is empty"<<endl;;
    hist->GetXaxis()->SetTitle("M_{T*} [GeV]");
    hist->GetYaxis()->SetTitle("events");
    hist->SetTitle(channel + ": " + histname);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(i+1);
    hist->SetLineColor(i+1);
    hist->Draw("hist same"); 
  }

  c1_hist->SaveAs(channel+"_M-"+masspoint+"_"+histname+".pdf");
}
