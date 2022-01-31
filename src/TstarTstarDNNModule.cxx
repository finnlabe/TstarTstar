#include <iostream>
#include <memory>
#include <string>

// UHH2 stuff
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/PhotonIds.h"
#include <UHH2/common/include/MuonIds.h>
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/MCWeight.h"

// TstarTstar
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarDNNHists.h"
#include "UHH2/TstarTstar/include/TstarTstarVariationHists.h"
#include "UHH2/TstarTstar/include/TstarTstarDNNInputHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarAllGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"
#include "UHH2/TstarTstar/include/TstarTstarGenMatch.h"
#include "UHH2/TstarTstar/include/NeuralNetworkModules.h"

// other
#include "UHH2/HOTVR/include/HOTVRIds.h"
#include "TGraphAsymmErrors.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

// x, alpha, n, sigma, mean
double crystalball_function(double x, double alpha, double n, double sigma, double mean) {
  // evaluate the crystal ball function
  if (sigma < 0.)     return 0.;
  double z = (x - mean)/sigma;
  if (alpha < 0) z = -z;
  double abs_alpha = std::abs(alpha);
  // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
  // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
  // double N = 1./(sigma*(C+D));
  if (z  > - abs_alpha)
    return std::exp(- 0.5 * z * z);
  else {
    //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
    double nDivAlpha = n/abs_alpha;
    double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
    double B = nDivAlpha -abs_alpha;
    double arg = nDivAlpha/(B-z);
    return AA * std::pow(arg,n);
  }
}

double bgest_polynom(double x, double p0, double p1, double p2) {
  return p0 + p1 * x + p2 * x * x;
}

/** \brief Module for the T*T*->ttbar gg MC based study
 *
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class TstarTstarDNNModule: public AnalysisModule {
public:

  explicit TstarTstarDNNModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  // ###### Modules ######
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<NeuralNetworkIncluder> DNN_Includer;


  // ##### Histograms #####
  std::unique_ptr<Hists> h_topcheck, h_topcheck_reweighted, h_AfterDNNcut_06_UGLYFIX;
  std::unique_ptr<Hists> h_STreweighted, h_crosscheck;

  std::unique_ptr<Hists> h_newTaggerSR, h_newTaggerCR, h_newTagger_btagCR;

  std::unique_ptr<Hists> h_AfterDNNcut_02, h_AfterDNNcut_03, h_AfterDNNcut_04, h_AfterDNNcut_05, h_AfterDNNcut_06, h_AfterDNNcut_07, h_AfterDNNcut_08;
  std::unique_ptr<Hists> h_notDNNcut_02,   h_notDNNcut_03,   h_notDNNcut_04,   h_notDNNcut_05,   h_notDNNcut_06,   h_notDNNcut_07,   h_notDNNcut_08;

  std::unique_ptr<Hists> h_BtagControl_AfterDNNcut_02, h_BtagControl_AfterDNNcut_03, h_BtagControl_AfterDNNcut_04, h_BtagControl_AfterDNNcut_05, h_BtagControl_AfterDNNcut_06, h_BtagControl_AfterDNNcut_07, h_BtagControl_AfterDNNcut_08;
  std::unique_ptr<Hists> h_BtagControl_notDNNcut_02,   h_BtagControl_notDNNcut_03,   h_BtagControl_notDNNcut_04,   h_BtagControl_notDNNcut_05,   h_BtagControl_notDNNcut_06,   h_BtagControl_notDNNcut_07,   h_BtagControl_notDNNcut_08;

  std::unique_ptr<Hists> h_DNN, h_DNN_newTagger, h_DNN_newTagger_2, h_DNN_BtagControl, h_DNN_reweighted, h_DNN_reweighted_2;
  std::unique_ptr<Hists> h_DNN_lowST, h_DNN_highST, h_DNN_lowDNN, h_DNN_highDNN, h_DNN_highST_lowDNN, h_DNN_highST_highDNN;
  std::unique_ptr<Hists> h_AfterDNN, h_AfterDNN_lowST, h_AfterDNN_highST, h_AfterDNN_lowDNN, h_AfterDNN_highDNN, h_AfterDNN_highST_lowDNN, h_AfterDNN_highST_highDNN;

  std::unique_ptr<Hists> h_highLepton, h_highLepton_AfterDNNcut_06, h_highLepton_notDNNcut_06;

  std::unique_ptr<Hists> h_SFVariations;

  // ###### Control switches ######
  bool debug = false;
  bool do_masspoint = false;
  bool doAddInputs = false;


  // ###### Handles ######
  uhh2::Event::Handle<int> h_flag_toptagevent;
  uhh2::Event::Handle<int> h_flag_muonevent;
  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<bool> h_is_btagevent;


  uhh2::Event::Handle<double> h_ST_weight;

  uhh2::Event::Handle<double> h_DNN_output;
  uhh2::Event::Handle<bool> h_do_masspoint;
  uhh2::Event::Handle<double> h_ST;
  uhh2::Event::Handle<bool> h_DoAddInputs;
  uhh2::Event::Handle<double> h_newTagger;
  uhh2::Event::Handle<TString> h_region;


  // ###### other parameters ######
  bool is_MC;
  bool is_TTbar = false;
  bool is_Signal = false;

  bool is_datadriven_BG_run = false;

  TGraphAsymmErrors* bgest_purity;

};


TstarTstarDNNModule::TstarTstarDNNModule(Context & ctx){

  // setting debug from xml file
  if(ctx.get("debug", "<not set>") == "true") debug = true;

  // debug message
  if(debug) {
    cout << "Hello World from TstarTstarDNNModule!" << endl;
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
  }

  // ###### 0. Setting variables ######
  // MC or real data
  is_MC = ctx.get("dataset_type") == "MC";

  if(!is_MC) is_datadriven_BG_run = ctx.get("use_data_for") == "background_extrapolation"; // get which running mode to use for data

  if(is_MC) {
    is_TTbar = (ctx.get("dataset_version").find("TT") != std::string::npos);
    is_Signal = (ctx.get("dataset_version").find("Tstar") != std::string::npos);
  }

  // ###### 1. set up modules ######
  // primlep
  reco_primlep.reset(new PrimaryLepton(ctx));

  // DNN modules
  DNN_Includer.reset(new NeuralNetworkIncluder(ctx, do_masspoint));


  // 3. Set up Hists classes:
  h_topcheck.reset(new TstarTstarAllGenHists(ctx, "topcheck"));
  h_topcheck_reweighted.reset(new TstarTstarAllGenHists(ctx, "topcheck_reweighted"));

  h_crosscheck.reset(new TstarTstarHists(ctx, "crosscheck"));
  h_STreweighted.reset(new TstarTstarHists(ctx, "STreweighted"));

  h_newTaggerSR.reset(new TstarTstarHists(ctx, "newTaggerSR"));
  h_newTaggerCR.reset(new TstarTstarHists(ctx, "newTaggerCR"));
  h_newTagger_btagCR.reset(new TstarTstarHists(ctx, "newTagger_btagCR"));

  h_SFVariations.reset(new TstarTstarVariationHists(ctx, "SFVariations"));

  /**
  h_AfterDNNcut_02.reset(new TstarTstarHists(ctx, "AfterDNNcut_02"));
  h_notDNNcut_02.reset(new TstarTstarHists(ctx, "notDNNcut_02"));
  h_AfterDNNcut_03.reset(new TstarTstarHists(ctx, "AfterDNNcut_03"));
  h_notDNNcut_03.reset(new TstarTstarHists(ctx, "notDNNcut_03"));
  h_AfterDNNcut_04.reset(new TstarTstarHists(ctx, "AfterDNNcut_04"));
  h_notDNNcut_04.reset(new TstarTstarHists(ctx, "notDNNcut_04"));
  h_AfterDNNcut_05.reset(new TstarTstarHists(ctx, "AfterDNNcut_05"));
  h_notDNNcut_05.reset(new TstarTstarHists(ctx, "notDNNcut_05"));
  h_AfterDNNcut_06.reset(new TstarTstarHists(ctx, "AfterDNNcut_06"));
  h_notDNNcut_06.reset(new TstarTstarHists(ctx, "notDNNcut_06"));
  h_AfterDNNcut_07.reset(new TstarTstarHists(ctx, "AfterDNNcut_07"));
  h_notDNNcut_07.reset(new TstarTstarHists(ctx, "notDNNcut_07"));
  h_AfterDNNcut_08.reset(new TstarTstarHists(ctx, "AfterDNNcut_08"));
  h_notDNNcut_08.reset(new TstarTstarHists(ctx, "notDNNcut_08"));
  **/

  h_AfterDNNcut_06_UGLYFIX.reset(new TstarTstarHists(ctx, "AfterDNNcut_06_UGLYFIX"));

  /**
  h_BtagControl_AfterDNNcut_02.reset(new TstarTstarHists(ctx, "AfterDNNcut_02_BtagControl"));
  h_BtagControl_notDNNcut_02.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_02"));
  h_BtagControl_AfterDNNcut_03.reset(new TstarTstarHists(ctx, "BtagControl_AfterDNNcut_03"));
  h_BtagControl_notDNNcut_03.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_03"));
  h_BtagControl_AfterDNNcut_04.reset(new TstarTstarHists(ctx, "BtagControl_AfterDNNcut_04"));
  h_BtagControl_notDNNcut_04.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_04"));
  h_BtagControl_AfterDNNcut_05.reset(new TstarTstarHists(ctx, "BtagControl_AfterDNNcut_05"));
  h_BtagControl_notDNNcut_05.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_05"));
  h_BtagControl_AfterDNNcut_06.reset(new TstarTstarHists(ctx, "BtagControl_AfterDNNcut_06"));
  h_BtagControl_notDNNcut_06.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_06"));
  h_BtagControl_AfterDNNcut_07.reset(new TstarTstarHists(ctx, "BtagControl_AfterDNNcut_07"));
  h_BtagControl_notDNNcut_07.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_07"));
  h_BtagControl_AfterDNNcut_08.reset(new TstarTstarHists(ctx, "BtagControl_AfterDNNcut_08"));
  h_BtagControl_notDNNcut_08.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_08"));
  **/

  h_highLepton.reset(new TstarTstarHists(ctx, "highLepton"));
  h_highLepton_AfterDNNcut_06.reset(new TstarTstarHists(ctx, "highLepton_AfterDNNcut_06"));
  h_highLepton_notDNNcut_06.reset(new TstarTstarHists(ctx, "highLepton_notDNNcut_06"));

  h_DNN.reset(new TstarTstarDNNHists(ctx, "DNN"));
  h_DNN_newTagger.reset(new TstarTstarDNNHists(ctx, "DNN_newTagger"));
  h_DNN_BtagControl.reset(new TstarTstarDNNHists(ctx, "DNN_BtagControl"));
  h_DNN_reweighted.reset(new TstarTstarDNNHists(ctx, "DNN_reweighted"));

  h_AfterDNN.reset(new TstarTstarHists(ctx, "AfterDNN"));

  // ###### 4. init handles ######
  h_evt_weight = ctx.get_handle<double>("evt_weight");
  h_flag_toptagevent = ctx.get_handle<int>("flag_toptagevent");
  h_DoAddInputs = ctx.declare_event_output<bool>("doAddInputs");
  h_newTagger = ctx.declare_event_output<double>("newTagger");
  h_is_btagevent = ctx.get_handle<bool>("is_btagevent");
  h_ST = ctx.get_handle<double>("ST");
  h_region = ctx.declare_event_output<TString>("region");

  h_ST_weight = ctx.declare_event_output<double>("ST_weight");

  h_DNN_output = ctx.get_handle<double>("DNN_output");

  TFile *f = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/CMSSW_10_2_17/src/UHH2/TstarTstar/macros/files/bgest_purity.root");
  bgest_purity = (TGraphAsymmErrors*)f->Get("purity");

}


bool TstarTstarDNNModule::process(Event & event) {

  // debug message
  if(debug){cout << endl << "TstarTstarDNNModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  // fixing poitential empty handle
  event.set(h_region,"not set");

  // reapply weights
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  // get ST weights
  double ST_weight = event.get(h_ST_weight);

  // set primary lepton
  reco_primlep->process(event);

  // setting addInputs
  event.set(h_DoAddInputs, doAddInputs);

  // placeholder as this flag seems to be wrong!!!
  int N_jets_btag_loose = 0;
  int N_jets_btag_medium = 0;
  int N_jets_btag_tight = 0;
  bool pass_btagcut = false;
  for(const auto & jet : *event.jets) {
    if(jet.btag_DeepCSV() > 0.2219) N_jets_btag_loose++;
    if(jet.btag_DeepCSV() > 0.2219) pass_btagcut = true;
    if(jet.btag_DeepCSV() > 0.6324) N_jets_btag_medium++;
    if(jet.btag_DeepCSV() > 0.8958) N_jets_btag_tight++;
  }

  h_crosscheck->fill(event);
  event.weight = ST_weight * event.get(h_evt_weight);
  if(pass_btagcut) h_STreweighted->fill(event);
  event.weight = event.get(h_evt_weight);


  // ################
  // ### DNN Part ###
  // ################

  // including trained DNN model
  if(debug) cout << "Include DNN model" << endl;
  DNN_Includer->process(event);

  if(debug) cout << "Done DNN include" << endl;
  if(event.get(h_DNN_output) > 0.6) h_AfterDNNcut_06_UGLYFIX->fill(event);

  if(debug) cout << "Start filling hists" << endl;
  // hists
  if(pass_btagcut) {
    h_DNN->fill(event);
    h_topcheck->fill(event);
    if(is_TTbar) event.weight *= ST_weight;
    h_DNN_reweighted->fill(event);
    h_topcheck_reweighted->fill(event);
    event.weight = event.get(h_evt_weight);

    // some more plotting
    double DNNoutput = event.get(h_DNN_output);
    h_AfterDNN->fill(event);

  }

  /**
  else {
    h_DNN_BtagControl->fill(event);
    if(event.get(h_DNN_output) > 0.2) h_BtagControl_AfterDNNcut_02->fill(event);
    else h_BtagControl_notDNNcut_02->fill(event);
    if(event.get(h_DNN_output) > 0.3) h_BtagControl_AfterDNNcut_03->fill(event);
    else h_BtagControl_notDNNcut_03->fill(event);
    if(event.get(h_DNN_output) > 0.4) h_BtagControl_AfterDNNcut_04->fill(event);
    else h_BtagControl_notDNNcut_04->fill(event);
    if(event.get(h_DNN_output) > 0.5) h_BtagControl_AfterDNNcut_05->fill(event);
    else h_BtagControl_notDNNcut_05->fill(event);
    if(event.get(h_DNN_output) > 0.6) h_BtagControl_AfterDNNcut_06->fill(event);
    else h_BtagControl_notDNNcut_06->fill(event);
    if(event.get(h_DNN_output) > 0.7) h_BtagControl_AfterDNNcut_07->fill(event);
    else h_BtagControl_notDNNcut_07->fill(event);
    if(event.get(h_DNN_output) > 0.8) h_BtagControl_AfterDNNcut_08->fill(event);
    else h_BtagControl_notDNNcut_08->fill(event);
  }
  **/

  if(debug) cout << "Defining CB function 1" << endl;

  // Additional decorrelation through "varying cut" on DNN output
  // the function used here was obtained by a fit in the macro "decorrelatedTagger.C"
  // 1  Constant     4.65951e-01   1.36602e-03  -9.94180e-06  -2.24498e-01
  // 2  Mean         7.13578e+02   4.00979e+00   1.37526e-03   6.80281e-05
  // 3  Sigma        2.40880e+02   7.13853e+00   8.82481e-03  -1.99541e-04
  // 4  Alpha       -1.43648e-01   6.23610e-03   1.21231e-05  -3.65818e-02
  // 5  N            5.14847e+05   2.62580e+05  -5.73867e+01   1.95656e-10
  double crystal_constant = 8.09337e-01;
  double crystal_mean = 6.86757e+02;
  double crystal_sigma = 4.37874e+02;
  double crystal_alpha = -2.14474e-01;
  double crystal_n = 7.97461e+05;
  double secondPart = 1 - (crystal_constant * crystalball_function(event.get(h_ST), crystal_alpha, crystal_n, crystal_sigma, crystal_mean));
  double newTagger = event.get(h_DNN_output) - secondPart;
  event.set(h_newTagger, newTagger);

  // datadriven background estimation
  if(is_datadriven_BG_run) {
    if(debug) cout << "Doing datadriven BG estimation" << endl;
    if(pass_btagcut) {
      pass_btagcut = false; // in this running scheme, we won't use the data that actually would go into the SR or CR
    } else {
      pass_btagcut = true; // we will use this data however!
      newTagger = -1; // for the moment, use all data for the CR

      // definition of the transfer function
      // NAME      VALUE            ERROR          SIZE      DERIVATIVE
      // p0       -1.20025e+00   3.82700e-01   1.31302e-04   2.81236e-04
      // p1        4.71172e-03   5.79365e-04   1.02640e-07   4.94329e-01
      // p2       -3.70059e-07   2.01813e-07   6.21146e-11   1.00150e+03
      double p0 = -6.71252e-01;
      double p1 = 4.16795e-03;
      double p2 = 1.86044e-07;
      if(true) cout << "ST: " << event.get(h_ST) << endl;
      double transfer_value = bgest_polynom(event.get(h_ST), p0, p1, p2);
      if(true) cout << "transfer value: " << transfer_value << endl;
      double purity_value = bgest_purity->Eval(event.get(h_ST));
      if(true) cout << "purity: " << purity_value << endl;
      event.weight *= transfer_value*purity_value;
    }
  }

  if(pass_btagcut) {
    h_DNN_newTagger->fill(event);
    if(newTagger > 0) {
      h_newTaggerSR->fill(event);
      event.set(h_region, "SR");
      h_SFVariations->fill(event);
    }
    else {
      h_newTaggerCR->fill(event);
      event.set(h_region,"CR1");
    }
  }
  else {
    if(newTagger > 0) {
      h_newTagger_btagCR->fill(event);
      event.set(h_region,"CR2");
    }
  }

  // testing tighter pt cut
  if(debug) cout << "testing tighter pt cut" << endl;
  if(pass_btagcut) {
    if(event.muons->size() == 1) {
      h_highLepton->fill(event);
      if(event.get(h_DNN_output) > 0.6) h_highLepton_AfterDNNcut_06->fill(event);
      else h_highLepton_notDNNcut_06->fill(event);
    }
    else {
      if(event.electrons->at(0).pt() > 50 && event.met->pt() > 80) {
        h_highLepton->fill(event);
        if(event.get(h_DNN_output) > 0.6) h_highLepton_AfterDNNcut_06->fill(event);
        else h_highLepton_notDNNcut_06->fill(event);
      }
    }

  }

  if(debug){cout << "Done ##################################" << endl;}
  return true;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarDNNModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarDNNModule)

}
