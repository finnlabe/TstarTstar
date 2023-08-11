

void backgroundEstimation_compareFunctions () {

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/";
  TString filename_base = "alphaFunction";

  std::vector<TString> variations = {"JEC_down", "JEC_up"};


  TCanvas *canvas = new TCanvas("canvas", "c", 600, 500);

  TFile *input_baseline = TFile::Open(path+filename_base+".root");
  TF1 *function_baseline = (TF1*)((TF1*) input_baseline->Get("fit2"))->Clone("baseline");

  function_baseline->Draw();

  for (auto variation : variations) {
    TFile *input = TFile::Open(path+filename_base+"_"+variation+".root");
    TF1 *function = (TF1*)((TF1*) input_baseline->Get("fit2"))->Clone(variation);

    TF1 ratio = TF1(variation+"_ratio", variation+"/baseline", 0, 4000);

    function->SetLineColor(2);
    function->Draw("same");

  }

  canvas->SaveAs("plots/BGest_variation.pdf");

}
