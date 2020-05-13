// Efficiency histograms for trigger related studies
// author: F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q efficiencies.C

void efficiencies(){
    
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
  Double_t w = 800;
  Double_t h = 600;
  
  // Reuseable buffers
  TFile *input;
  TH1D *hist;

  // Drawing Definitions
  TCanvas *canvas = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.22,0.175,0.5,0.33);
  //canvas->SetLogy();

  // Defining paths
  std::vector<TString> paths = {"/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Selection/hadded/","/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/MCStudy/hadded/"};
  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";
  TString histname = "N_jets";

  // Defining Steps
  std::vector<std::vector<TString>> steps = {{"After2D"}, {"AfterReco_Full"}};

  // Defining Samples
  std::vector<TString> samples = {"TstarTstar_M-700.root","TstarTstar_M-800.root", "TstarTstar_M-900.root", "TstarTstar_M-1000.root", "TstarTstar_M-1100.root", "TstarTstar_M-1200.root", "TstarTstar_M-1300.root", "TstarTstar_M-1400.root", "TstarTstar_M-1500.root", "TstarTstar_M-1600.root"};

  // Defining Drawing options
  std::vector<TString> labels = {"Preselection", "Reconstruction"};

  // ########################
  // ## Finish Definitions ##
  // ########################

  // Saving initial values
  TString initial_step = "AfterCommon";
  TString initial_path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Preselection/hadded/";
  std::vector<double> initial;
  for(const auto & sample : samples){
    input = TFile::Open(initial_path+fileprefix+sample);
    hist = (TH1D*)input->Get(initial_step+"/"+histname);
    initial.push_back(hist->Integral());
  }

  // Saving Vector for step values
  std::vector<std::vector<double>> values;
  int count_path = 0;
  for(const auto & path : paths){
    for(const auto & step : steps.at(count_path)){
      std::vector<double> values_sample;
      for(const auto & sample : samples){
	input = TFile::Open(path+fileprefix+sample);
	hist = (TH1D*)input->Get(step+"/"+histname);
	values_sample.push_back(hist->Integral());
      }
      values.push_back(values_sample);
    }
    count_path++;
  }

  for(const auto & valvec : values){
    std::cout << "Step" << endl;
    int i = 0;
    for(const auto & value : valvec){
      std::cout << value/initial.at(i) << " " << endl;
      i++;
    }
  }

}
