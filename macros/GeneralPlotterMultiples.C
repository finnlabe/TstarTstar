// Just a basic macro to create simple plots from given hists.
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q GeneralPlotterMultiples.C

void GeneralPlotterMultiples(){

  TString path_pre = "/nfs/dust/cms/user/flabe/TstarTstar/data/";
  std::vector<TString> files = {"GenInfo/hadded/uhh2.AnalysisModuleRunner.MC.TstarTstar_M-700.root", "GenInfo/hadded/uhh2.AnalysisModuleRunner.MC.TstarTstar_M-1300.root", "GenInfo/hadded/uhh2.AnalysisModuleRunner.MC.TstarTstar_M-2000.root"};
  std::vector<TString> steps = {"GenRecoMatchedHists", "GenRecoMatchedHists", "GenRecoMatchedHists"};
  std::vector<TString> hists = {"topmatches", "topmatches", "topmatches"};
  std::vector<TString> labels = {"T* M-700", "T* M-1300", "T* M-2000"};
  bool normalize = true;

  Double_t w = 800;
  Double_t h = 600;

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  auto legend = new TLegend(0.1,0.8,0.3,0.9);

  for(int i = 0; i < files.size(); i++){
    std::cout << "Plotting: " << files.at(i) << " // " << steps.at(i) << " // " << hists.at(i) << endl;
    TFile *input = TFile::Open(path_pre+files.at(i));
    TH1D *hist = (TH1D*)input->Get(steps.at(i)+"/"+hists.at(i)); //histogram
    if(!hist) cout<<"Hist is empty"<<endl;;
    if(normalize) hist->Scale(1/hist->Integral());
    hist->SetMarkerStyle(0);
    hist->SetMarkerColor(i+1);
    hist->SetLineStyle(i+1);
    hist->GetYaxis()->SetRangeUser(0, .8);
    hist->Draw("hist same");
    legend->AddEntry(hist, labels.at(i));
  }

  legend->Draw("same");

  c1_hist->SaveAs("plot.pdf");
}
