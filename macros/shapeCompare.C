// Just a basic macro to create simple plots from given hists.
// author: F. Labe
// Run it with following command:
// root -l -b -q shapeCompare.C

void shapeCompare(){

  gStyle->SetOptStat(0);

  Double_t w = 800;
  Double_t h = 600;

  TString histname = "pt_ST";
  std::vector<TString> files = {
    "/nfs/dust/cms/user/flabe/TstarTstar/data/Analysis/hadded/uhh2.AnalysisModuleRunner.MC.TTbar.root",
    "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/hadded/uhh2.AnalysisModuleRunner.MC.TTbar.root",
    "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/hadded/uhh2.AnalysisModuleRunner.MC.TTbar.root",
    "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/hadded/uhh2.AnalysisModuleRunner.MC.TTbar.root",
  };
  std::vector<TString> folders = {"main", "AfterDNNcut_02", "AfterDNNcut_05", "AfterDNNcut_08"};
  assert(file.size() == folder.size());

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.22,0.175,0.6,0.40);

  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  c1_hist->SetLogy();

  int i = 0;
  for (auto file : files) {
    cout << "Processing " << file << "." << endl;
    TString folder = folders.at(i);
    TFile *input = TFile::Open(file);
    if(!input) cout<<"File is empty"<<endl;
    TH1D *hist = (TH1D*)input->Get(folder+"/"+histname); //histogram
    if(!hist) cout<<"Hist is empty"<<endl;
    hist->Scale(1/hist->Integral());
    hist->GetYaxis()->SetTitle("events");
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(i+1);
    hist->SetLineColor(i+1);
    if(i == 0) hist->Draw("hist");
    else hist->Draw("hist same");
    leg->AddEntry(hist, folder, "l");
    i++;
  }

  leg->Draw("same");
  c1_hist->SaveAs("plots/shape_"+histname+".pdf");
}
