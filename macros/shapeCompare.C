// Just a basic macro to create simple plots from given hists.
// author: F. Labe
// Run it with following command:
// root -l -b -q shapeCompare.C

void shapeCompare(){

  Double_t w = 800;
  Double_t h = 600;

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/";
  TString filename = "FullSignal.root";
  TString histname = "pt_ST_jets";

  TString subpathA = "/Analysis/";
  TString folderA = "beforeReco";
  TString subpathB = "/DNN/";
  TString folderB = "AfterDNNcut_06"; // sure?

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  c1_hist->SetLogy();

  TFile *input1 = TFile::Open(path+subpathA+"/"+filename);
  if(!input1) cout<<"File A is empty"<<endl;
  TH1D *hist1 = (TH1D*)input1->Get(folderA+"/"+histname); //histogram
  if(!hist1) cout<<"Hist A is empty"<<endl;
  hist1->Scale(1/hist1->Integral());
  hist1->GetXaxis()->SetTitle("S_{T} [GeV]");
  hist1->GetYaxis()->SetTitle("events");
  hist1->SetMarkerStyle(20);
  hist1->SetMarkerColor(1);
  hist1->SetLineColor(1);
  hist1->Draw("hist");

  TFile *input2 = TFile::Open(path+subpathB+"/"+filename);
  TH1D *hist2 = (TH1D*)input2->Get(folderB+"/"+histname); //histogram
  if(!hist2) cout<<"Hist A is empty"<<endl;
  hist2->Scale(1/hist2->Integral());
  hist2->GetXaxis()->SetTitle("S_{T} [GeV]");
  hist2->GetYaxis()->SetTitle("events");
  hist2->SetMarkerStyle(20);
  hist2->SetMarkerColor(1);
  hist2->SetLineColor(2);
  hist2->Draw("hist same");

  c1_hist->SaveAs("shapeST_sig.pdf");
}
